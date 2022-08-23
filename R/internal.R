# Constructor to allow logging in with username/password
#' Instantiate the MGnifyR client object
#'
#' All functions in the MGnifyR package take a \code{mgnify_client} object as their first argument. While not essential
#' to querying the raw MGnify API (which is exposed as relative standard JSONAPI), the object allows the simple handling of both
#' user authentication and access to private data, and local on-disk caching of results.
#'
#' @importFrom httr POST
#' @importFrom httr content
#'
#' @param url (To be described)
#' @param username optional username to authenticate.
#' @param password optional password for authentication.
#' @param usecache whether to enable on-disk caching of results during this session. In most use cases should be TRUE.
#' @param cache_dir specifies a folder to contain the local cache. If NULL, and usecache is TRUE, the new subdirectory \code{.MGnifyR_cache}
#' in the current working directory will be used. Note that cached files are persistent, so the cache directory may be reused between sessions,
#' taking advantage of previously downloaded results. The directory will be created if it doesn't exist already.
#' @param warnings debug flag to print extra output during invocation of some MGnifyR functions. Defaults to FALSE.
#' @param use_memcache flag to indicate whether functional results obtained when \code{bulk_dl} is \code{TRUE} in \code{mgnify_get_analyses_results} should
#' be stored in an in-memory cache, rather than the cached input being re-read for each accession. this is currently NOT working
#' properly and should therefore be set \code{FALSE} (the default). It has the potential to speed up searches considerably though, especially
#' for studies with a large number of samples, so will be implemented properly in the future.
#' @examples
#' my_client <- mgnify_client(username="Webin-1122334", password="SecretPassword", usecache=T, cache_dir = "/scratch/MGnify_cache_location")
#' @export
mgnify_client <- function(url=NULL,username=NULL,password=NULL,usecache=F,cache_dir=NULL, warnings=F, use_memcache=F){
    if (is.null(url)){
        url <- baseurl
    }

    authtok <- NA_character_

    #Check to see if we're goint to try and get an authentication token:
    if (!is.null(username) && !is.null(password)){
        r <- POST(paste(url, "utils/token/obtain", sep="/"),
                                     body=list(username=username, password=password),
                                     encode="json")
        cont <- content(r)
        if ("data" %in% names(cont)){
            authtok <- cont$data$token
        }
        else{
            stop("Failed to authenticate")
        }
    }
    #Assume we're not using it
    cachepath <- NA_character_
    if(usecache){
        if (is.null(cache_dir) ){
            cachepath <- paste(getwd(),'.MGnifyR_cache',sep="/")
        }else{
            cachepath <- cache_dir
        }
        #Make it if needed - assume the user is sensible and the path will work...
        dir.create(cachepath,showWarnings = F)
    }

    #Return the final object
    #@importFrom methods new
    new("mgnify_client", url=url, authtok=authtok, cache_dir = cachepath, warnings=warnings, memcache=list(), use_memcache=use_memcache)
}

#Retrieves combined study/sample/analysis metadata - not exported
mgnify_get_single_analysis_metadata <- function(client=NULL, accession, usecache=T, maxhits=-1){

    dat <- mgnify_retrieve_json(client, paste("analyses", accession, sep="/"), usecache = usecache, maxhits = maxhits)
    #There ~should~ be just a single result
    top_data <- dat[[1]]
    analysis_df <- mgnify_attr_list_to_df_row(top_data, metadata_key = "analysis-summary")

    #cat(str(analysis_metadata))
    #Build up the metadata dataframe from the analyses_metadata_headers vector:
    sample_met <- mgnify_retrieve_json(client, complete_url = top_data$relationships$sample$links$related, usecache = usecache)
    study_met <- mgnify_retrieve_json(client, complete_url = top_data$relationships$study$links$related, usecache = usecache)

    sample_df <- mgnify_attr_list_to_df_row(sample_met[[1]], metadata_key = "sample-metadata")
    #It turns out that a sample might not be part of a study - if it's been harvested...
    #So tryCatch it and return an empy df row if things go south.
    study_df <- tryCatch(mgnify_attr_list_to_df_row(study_met[[1]]), error=function(X) {
        warning(paste("Failed to find study metadata for ", accession, sep=""))
        data.frame(accession=NA)
    })

    colnames(sample_df) <- paste("sample",colnames(sample_df), sep="_")
    colnames(study_df) <- paste("study",colnames(study_df), sep="_")
    colnames(analysis_df) <- paste("analysis",colnames(analysis_df), sep="_")

    rownames(sample_df) <- rownames(analysis_df)
    rownames(study_df) <- rownames(analysis_df)
    full_df <- cbind(analysis_df, study_df, sample_df)

    #extras - include some more metadata from various places
    #assembly accesion
    if("id" %in% names(top_data$relationships$assembly$data)){
        full_df$assembly_accession <- top_data$relationships$assembly$data$id
    }
    #run accession
    if("id" %in% names(top_data$relationships$run$data)){
        full_df$run_accession <- top_data$relationships$run$data$id
    }

    #biom (from the sample metadata)
    tryCatch({
        full_df$biome_string <- sample_met[[1]]$relationships$biome$data$id
    },
        error <- function(x) warning("Error finding biome entry")
    )

    full_df
}

#' @importFrom urltools parameters
#' @importFrom httr GET
#' @importFrom httr write_disk
#' @importFrom phyloseq import_biom
#' @importFrom phyloseq tax_table
#' @importFrom phyloseq parse_taxonomy_greengenes
#' @importFrom phyloseq sample_names
#' @importFrom phyloseq sample_data
#' @importFrom phyloseq phy_tree
#' @importFrom ape read.tree
#' @importFrom utils tail

#UPDATE ME SO THAT TREES (if available) GET GRABBED AS WELL!!!
# Not exported - get a single biom file and convert it to a phyloseq object.
mgnify_get_single_analysis_phyloseq <- function(client=NULL, accession, usecache=T, downloadDIR=NULL, tax_SU="SSU", get_tree=FALSE){
    metadata_df <- mgnify_get_single_analysis_metadata(client, accession, usecache=usecache)

    analysis_data <- mgnify_retrieve_json(client, paste("analyses",accession,sep="/"), usecache = usecache)
    download_url <- analysis_data[[1]]$relationships$downloads$links$related
    analysis_downloads <- mgnify_retrieve_json(client, complete_url = download_url,usecache = usecache)

    #Depending on the pipeline version, there may be more than one OTU table available (LSU/SSU), so try and get the
    #one specified in tax_SU - otherwise spit out a warning and grab the generic (older pipelines)
    available_biom_files <- analysis_downloads[grepl('JSON Biom', sapply(analysis_downloads, function(x) x$attributes$`file-format`$name))]
    biom_position <- grepl(tax_SU, sapply(available_biom_files, function(x) x$attributes$`group-type`))
    if(sum(biom_position) == 0){
        if(client@warnings){
        warning("Unable to locate requested taxonomy type ",tax_SU,". This is likely due to the current analysis having been performed on an older version of the MGnify pipeline.
                         The available BIOM file will be used instead.")
        }
        biom_url <- available_biom_files[[1]]$links$self
    }else{
        biom_url <- available_biom_files[biom_position][[1]]$links$self
    }

    #Can specify a seperate dir for saving biom files, otherwise they end up in the client@cachdir folder, under "bioms"
    if (is.null(downloadDIR)){
        downloadDIR <- paste(client@cache_dir,"biom_files",sep="/")
        dir.create(downloadDIR, recursive = T, showWarnings = client@warnings)
    }
    #Clear out any ?params after the main location - don't need them for this
    parameters(biom_url) <- NULL

    fname <- tail(strsplit(biom_url, '/')[[1]], n=1)
    biom_path <- paste(downloadDIR, fname, sep="/")

    ## Quick check to see if we should clear the disk cache ~for this specific call~ - used for debugging
    # and when MGnify breaks
    if(usecache & client@clear_cache){
        message(paste("clear_cache is TRUE: deleting ",biom_path, sep=""))
        tryCatch(unlink(biom_path), error=warning)
    }

    if (! file.exists(biom_path)){#} | !use_downloads ){
        GET(biom_url, write_disk(biom_path, overwrite = T))
    }
    #Load in the phlyloseq object
    psobj <- import_biom(biom_path)
    #Need to check if the taxonomy was parsed correctly - depending on the pipeline it may need a bit of help:
    if (ncol(tax_table(psobj)) == 1){
        psobj <- import_biom(biom_path, parseFunction = parse_taxonomy_qiime)
    }
    if(! "Kingdom" %in% colnames(tax_table(psobj))){
        psobj <- import_biom(biom_path, parseFunction = parse_taxonomy_greengenes)
    }

    #The biom files have a single column of unknown name - maybe it's the original sample name?.
    # It's rewritten as sample_run_analysis accession, with
    # the original value stored in the sample_data (just in case it's needed later)
    orig_samp_name <- sample_names(psobj)[[1]]
    newsampname <- rownames(metadata_df)[1]
    metadata_df[1,"orig_samp_name"] <- orig_samp_name
    sample_names(psobj) <- newsampname
    sample_data(psobj) <- metadata_df

    #Finally, do we want to the phylogenetic tree? If so, is it there?
    if(get_tree){
        #is there a tree?
        tvec <- grepl('Phylogenetic tree', sapply(analysis_downloads, function(x) x$attributes$`description`$label))
        if(any(tvec)){
            tree_url <- analysis_downloads[tvec][[1]]$links$self
            #Clear out any ?params after the main location - don't need them for this
            parameters(tree_url) <- NULL
            
            fname <- tail(strsplit(tree_url, '/')[[1]], n=1)
            tree_path <- paste(downloadDIR, fname, sep="/")

            ## Quick check to see if we should clear the disk cache ~for this specific call~ - used for debugging
            # and when MGnify breaks
            if(usecache & client@clear_cache){
                message(paste("clear_cache is TRUE: deleting ",tree_path, sep=""))
                tryCatch(unlink(tree_path), error=warning)
            }

            if (! file.exists(tree_path)){#} | !use_downloads ){
                GET(tree_url, write_disk(tree_path, overwrite = T ))
            }
        }
        phylo_tree <- read.tree(tree_path)
        phy_tree(psobj) <- phylo_tree
    }
    psobj
}

#' @importFrom urltools parameters
#' @importFrom httr GET
#' @importFrom httr write_disk
#' @importFrom dplyr bind_rows
#' @importFrom utils tail
#' @importFrom utils read.csv2

#' types can be found using \code{ names(MGnifyR::analyses_results_type_parsers)}. Note that not depending on the particular analysis type, pipeline

#Retrieves combined study/sample/analysis metadata - not exported
mgnify_get_single_analysis_results <- function(client=NULL, accession, retrievelist=c(), usecache=T, maxhits=-1, bulk_files=F){
    metadata_df <- mgnify_get_single_analysis_metadata(client, accession, usecache=usecache, maxhits = maxhits)

    #Should we try and grab the study's full TSV download rather than parse through the JSON API? Doing
    #so has the potential to use a LOT more disk space, along with potentially increased data download. It should
    #be faster though, except in pathological cases (e.g. only 1 sample per 1000 sample study required). As with everything
    #else, we make use of local caching to speed things along.
    if(bulk_files){
        downloadDIR <- paste(client@cache_dir, "tsv", sep="/")
        if(!dir.exists(downloadDIR)){
            dir.create(downloadDIR, recursive = T, showWarnings = client@warnings)
        }

        available_downloads_json <- mgnify_retrieve_json(client, path=paste("studies",metadata_df$study_accession,"downloads", sep="/"), usecache = usecache)
        #cat(str(dput(available_downloads_json)))
        parsed_results <- lapply(available_downloads_json, function(r) {
            #Figure out the label mapping
            cur_lab <- r$attributes$description$label

            cur_pipeversion <- r$relationships$pipeline$data$id

            #There MUST be a better way to do the line below...
            cur_type <- names(analyses_results_bulk_file_names[is.finite(match(analyses_results_bulk_file_names, cur_lab))])

            if (!identical(cur_type, character(0))){

                #Check the pipeline versions match
                if (cur_pipeversion == metadata_df$`analysis_pipeline-version`){
                    if( cur_type %in% retrievelist) {

                    #Get the url
                    data_url <- r$links$self

                    #Clear off extraneous gubbins
                    parameters(data_url) <- NULL

                    #build the cache filename
                    fname <- tail(strsplit(data_url, '/')[[1]], n=1)

                    #At this point we might have alread got the data we want loaded. Check the memory cache object

                    if(client@use_memcache & (cur_type %in% names(mgnify_memory_cache)) & (mgnify_memory_cache[cur_type]["fname"] == fname)){
                        tmp_df <- mgnify_memory_cache[cur_type][["data"]]
                    }else{
                        #Nope - gonna have to load it up from disk or grab it from t'interweb
                        data_path <- paste(downloadDIR, fname, sep="/")

                        if(usecache & client@clear_cache){
                            message(paste("clear_cache is TRUE: deleting ", data_path, sep=""))
                            tryCatch(unlink(data_path), error=warning)
                        }

                        if (! file.exists(data_path)){#} | !use_downloads ){
                            GET(data_url, write_disk(data_path, overwrite = T ))
                        }

                        #Load the file (might be big so save it in the 1 deep cache)
                        tmp_df <- read.csv2(data_path, sep="\t", header = T, stringsAsFactors = F)
                    }
                    #Save it in memory using "super assignment" - which I'm not really sure about but it seems to work... thing'd be
                    # much easier if R passed objects by reference.
                    if(client@use_memcache){
                        mgnify_memory_cache[[cur_type]] <<- list(data=tmp_df, fname=fname)
                    }

                    #Because there seem to be "mismatches" between the JSON and downloadable files, and maybe some issues with missing
                    #downloads, we have to check if we actually got a valid file:
                    if(ncol(tmp_df) <3){

                        warning(paste("Invalid download for", accession, sep=" "))
                        return(NULL)
                    }

                    #Need to figure out how many columns to keep - the first one is always an ID, but need to keep some others as well...
                    i <- 1
                    #tmp_df <- read.csv('~/.MGnify_cache/tsv/ERP108138_IPR_abundances_v4.1.tsv', sep="\t", header = T, stringsAsFactors = F)
                    while(any(is.na(suppressWarnings(as.numeric(tmp_df[,i] ))))){
                        i <- i+1
                    }
                    i <- i-1

                    #also need the column name for this particular analysis... As far as I can see they could be either assembly IDs or run ids. FFS.
                    #Assuming that both assembly and run won't be present...:
                    if("assembly_accession" %in% colnames(metadata_df)){
                        accession <- metadata_df$assembly_accession[[1]]
                    }else if("run_accession" %in% colnames(metadata_df)){
                        accession <- metadata_df$run_accession[[1]]
                    }
                #cat(accession)
                #    cat(str(head(tmp_df[1:5,1:5])))
                    #Break up the dataframe, only keeping the bits we need.

                    #so at this point we learn that some of the "download" files don't match the assembly/run IDs given in the JSON.
                    #For now, do a try/catch, chuck a warning, and then optionally go off and try again - this time from the JSON.
                    #No doubt this'll be fixed at some point in the future...

                    column_position <- match(accession, colnames(tmp_df))
                    if (is.na(column_position)){
                        warning(paste("Failed to find column",accession, sep=" "))
                        return(NULL)
                    }
                    keeper_columns <- c(seq(1,i), column_position)

                    #cat(keeper_columns)
                    tmp_df2 <- tmp_df[,keeper_columns]

                    tmp_colnames <- colnames(tmp_df2)
                    tmp_colnames[1] <- "accession"
                    tmp_colnames[length(tmp_colnames)] <- "count"

                    colnames(tmp_df2) <- tmp_colnames

                    tmp_df2$index_id <- tmp_df2$accession
                    rownames(tmp_df2) <- tmp_df2$accession
                    tmp_df2
                    }
                }
            }
        })
        #R is sometimes a bit ~awkward~
        names(parsed_results) <- names(analyses_results_bulk_file_names)[match(unlist(lapply(available_downloads_json,
                                                                                                                                                                                    function(x) x$attributes$description$label)),
                                                                                                                                                        analyses_results_bulk_file_names)]


    }else{
        #Now (re)load the analysis data:
        analysis_data <- mgnify_retrieve_json(client, paste("analyses",accession,sep="/"), usecache = usecache, maxhits = maxhits)
        #For now try and grab them all - just return the list - don't do any processing...
        all_results <- lapply(names(analyses_results_type_parsers), function(r) {
            if(r %in% retrievelist){
                tmp <- mgnify_retrieve_json(client, complete_url = analysis_data[[1]]$relationships[[r]]$links$related, usecache = usecache, maxhits=maxhits)
                tmp
            }
        })
        names(all_results) <- names(analyses_results_type_parsers)
        parsed_results <- sapply(names(all_results), function(x){
            all_json <- all_results[[x]]
            if(! is.null(all_json)){
                res_df <- do.call(bind_rows, lapply(all_json,analyses_results_type_parsers[[x]] ))
                rownames(res_df) <- res_df$index_id
                res_df
            }else{
                NULL
            }
        })
    }
    #Return the results...
    parsed_results
}

#' @importFrom mia loadFromBiom
#' @importFrom mia checkTaxonomy
#' @importFrom urltools parameters
#' @importFrom httr GET
#' @importFrom httr write_disk
#' @importFrom ape read.tree
#' @importFrom TreeSummarizedExperiment rowTree
#' @importFrom utils tail
#'
#' Get a single biom file and convert it to TreeSummarizedExperiment format
mgnify_get_single_analysis_treese <- function(client=NULL, accession, usecache=T, downloadDIR=NULL, tax_SU="SSU", get_tree=FALSE){

    metadata_df <- mgnify_get_single_analysis_metadata(client, accession, usecache=usecache)
    analysis_data <- mgnify_retrieve_json(client, paste("analyses",accession,sep="/"), usecache = usecache)
    download_url <- analysis_data[[1]]$relationships$downloads$links$related
    analysis_downloads <- mgnify_retrieve_json(client, complete_url = download_url,usecache = usecache)

    #Depending on the pipeline version, there may be more than one OTU table available (LSU/SSU), so try and get the
    #one specified in tax_SU - otherwise spit out a warning and grab the generic (older pipelines)
    available_biom_files <- analysis_downloads[grepl('JSON Biom', sapply(analysis_downloads, function(x) x$attributes$`file-format`$name))]
    biom_position <- grepl(tax_SU, sapply(available_biom_files, function(x) x$attributes$`group-type`))
    if(sum(biom_position) == 0){
        if(client@warnings){
            warning("Unable to locate requested taxonomy type ",tax_SU,". This is likely due to the current analysis having been performed on an older version of the MGnify pipeline. The available BIOM file will be used instead.")
        }
        biom_url <- available_biom_files[[1]]$links$self
    }else{
        biom_url <- available_biom_files[biom_position][[1]]$links$self
    }

    #Can specify a seperate dir for saving biom files, otherwise they end up in the client@cachdir folder, under "bioms"
    if (is.null(downloadDIR)){
        downloadDIR <- paste(client@cache_dir,"biom_files",sep="/")
        dir.create(downloadDIR, recursive = T, showWarnings = client@warnings)
    }
    #Clear out any ?params after the main location - don't need them for this
    parameters(biom_url) <- NULL

    fname <- tail(strsplit(biom_url, '/')[[1]], n=1)
    biom_path <- paste(downloadDIR, fname, sep="/")

    ## Quick check to see if we should clear the disk cache ~for this specific call~ - used for debugging
    # and when MGnify breaks
    if(usecache && client@clear_cache){
        message(paste("clear_cache is TRUE: deleting ",biom_path, sep=""))
        tryCatch(unlink(biom_path), error=warning)
    }

    if (! file.exists(biom_path)){#} | !use_downloads ){
        GET(biom_url, write_disk(biom_path, overwrite = T))
    }

    #Load in the TreeSummarizedExperiment object
    tse <- loadFromBiom(biom_path)

    #Need to check if the taxonomy was parsed correctly - depending on the pipeline it may need a bit of help:
    checkTaxonomy(tse)

    if(get_tree){
        #is there a tree?
        tvec <- grepl('Phylogenetic tree', sapply(analysis_downloads, function(x) x$attributes$`description`$label))
        if(any(tvec)){
            tree_url <- analysis_downloads[tvec][[1]]$links$self
            #Clear out any ?params after the main location - don't need them for this
            parameters(tree_url) <- NULL

            fname <- tail(strsplit(tree_url, '/')[[1]], n=1)
            tree_path <- paste(downloadDIR, fname, sep="/")

            ## Quick check to see if we should clear the disk cache ~for this specific call~ - used for debugging
            # and when MGnify breaks
            if(usecache && client@clear_cache){
                message(paste("clear_cache is TRUE: deleting ",tree_path, sep=""))
                tryCatch(unlink(tree_path), error=warning)
            }

            if (! file.exists(tree_path)){#} | !use_downloads ){
                GET(tree_url, write_disk(tree_path, overwrite = T ))
            }
        }

        row_tree <- read.tree(tree_path)
        rowTree(tse) <- row_tree
    }
    tse
}

# helper function for getting relative paths in the API
# Not everything is implemented here - just what we
# need to get to the download or run areas
# Given an accession x, we want to get the link to get the url for the
# corresponding typeY
#' \code{JSONAPI} path for child elements
#'
#' \code{mgnify_get_x_for_y} determines the location of \code{typeY} child objects of \code{x} (type \code{typeX})
#'
#' This helper function, principally intended to be used internally, is used to match up related objects within the path. The inherently
#' unhierarchical nature of the MGnify API makes it a bit inconsistent. This function acts as a quick way to determine how to get from
#' one type to another, without having to special case within the code.
#' @param client MGnifyR client API object
#' @param x Accession ID \code{char} of parent object
#' @param typeX Type of accession \code{x}
#' @param typeY Type of child object to return
#' @param usecache (To be described)
#' @return \code{char} complete url to access the result. Note this query is not run from here - just the URL is returned
#' @examples
#' cl <- new("mgnify_client")
#' mgnify_get_x_for_y(cl, "MGYS00005126", "studies", "samples")
## @export
mgnify_get_x_for_y <- function(client, x, typeX, typeY, usecache=F){
    #This one's easy - just rearrange the URLs
    #if(typeX=="samples" & typeY %in% c("runs","studies")){
    #    paste( typeX,x,typeY, sep="/")
    #}else if(typeX=="runs" & typeY == "analyses"){
    #    paste( typeX,x,typeY, sep="/")
    #}
    #else{
        #Do it the hard way with a callout
        json_dat <- mgnify_retrieve_json(client, paste(typeX, x, sep="/"), usecache = usecache)
        #cat(str(json_dat))
        #tgt_access = json_dat[[1]]$relationships[[typeY]]$data$id
        #tgt_type = json_dat[[1]]$relationships[[typeY]]$data$type
        #paste(tgt_type,tgt_access,sep="/")
        json_dat[[1]]$relationships[[typeY]]$links$related
        #substr(tgt_url, nchar(client@url) + 1, nchar(tgt_url))
    #}
}

# Internal function to actually perform the http request. Build up the URL then issues
# a GET, parsing the returned JSON into a nested list (uses \code{jsonlite} internally?)
# Previously cached results may be retrieved from disk without resorting to calling the MGnify server.

#' Low level MGnify API handler
#'
#' \code{mgnify_retrieve_json} deals with handles the actual HTTP GET calls for the MGnifyR package, handling API pagination,
#' local result caching, and authentication cookies for access
#' to restricted or pre-release datasets.Although principally intended for internal MGnifyR use , it's exported for direct invocation.
#' Generally though it's not recommended for use by users.
#'
#' @importFrom urltools parameters
#' @importFrom httr add_headers
#' @importFrom httr GET
#' @importFrom httr config
#' @importFrom httr content
#'
#' @param client MGnifyR client
#' @param path top level search point for the query. One of \code{biomes}, \code{samples}, \code{runs} etc. Basically includes
#' all parts of the URL between the base API url and the parameter specifications
#' @param complete_url \emph{complete} url to search, usuaally retrieved from previous query in the "related" section.
#' @param qopts named list or vector containing options/filters to be URL encoded and appended to query as key/value pairs
#' @param maxhits Maxmium number of data entries to return. The actual number of hits returned may be higher than this value,
#' as this parameter only clamps after each full page is processed. Set to <=0 to disable - i.e. retrieve all items.
#' @param usecache Should successful queries be cached on disk locally? There are unresolved questions about whether this is
#' a sensible thing to do, but it remains as an option. It probably makes sense for single accession grabs, but not for
#' (filtered) queries - which are liable to change as new data is added to MGnify. Also caching only works for the first page.
#' @param Debug Should we print out lots of information while doing the grabbing?
#' @return \code{list} of results after pagination is dealt with.
#' @export
mgnify_retrieve_json <- function(client, path="biomes", complete_url=NULL, qopts=NULL,maxhits=200, usecache = F, Debug=F){


    #client@warnings turns on debugging too:

    if(client@warnings){
        Debug <- T
    }
    # Set up the base url
    # Are we using internal paths?
    if (is.null(complete_url)){
        fullurl <- paste(client@url, path, sep="/")
    }
    #Or direct links from e.g. a "related" section
    else{
        #Set the full url, but clean off any existing parameters (page, format etc) as they'll be added back later:
        fullurl <- complete_url
        parameters(fullurl) <- NULL
        path <- substr(fullurl, nchar(client@url) + 2, nchar(fullurl))
    }

    #cat(fullurl)

    #convert to csv if filters are lists.
    #This doesn't check if they ~can~ be searched for in the API,
    #which is an issue since no error is returned by the JSON if the search
    #is invalid - we only get a result as if no query was present...
    tmpqopts <- lapply(qopts,function(x) paste(x,collapse = ','))

    #Include the json and page position options
    #full_qopts <- as.list(c(format="json", tmpqopts, page=1))
    full_qopts <- as.list(c(format="json", tmpqopts))
    #Build up the cache name anyway - even if it's not ultimately used:
    fname_list <- c(path, names(unlist(full_qopts)), unlist(full_qopts))
    cache_fname <- paste(fname_list,collapse = "_")
    cache_full_fname <- paste(client@cache_dir, '/', cache_fname, '.RDS', sep="")


    ## Quick check to see if we should clear the disk cache ~for this specific call~ - used for debugging
    # and when MGnify breaks
    if(usecache & client@clear_cache){
        message(paste("clear_cache is TRUE: deleting ", cache_full_fname, sep=""))
        tryCatch(unlink(cache_full_fname), error=warning)
    }

    # Do we want to try and use a cache to speed things up?
    if(usecache & file.exists(cache_full_fname)){
            final_data <- readRDS(cache_full_fname)
    }else{

        #Authorization: Bearer <your_token>
        if(!is.null(client@authtok)){
            add_headers(.headers = c(Authorization = paste("Bearer", client@authtok, sep=" ")))
        }
        res <- GET(url=fullurl, config(verbose=Debug), query=full_qopts )
        data <- content(res)

        #At this point, data$data is either a list of lists or a single named list. If it's a single entry, it needs embedding in
        #a list for consistency downstream
        #datlist is built up as a list of pages, where each entry must be another list. Thus, on the first page,
        #
        datlist <- list()
        if (!is.null(names(data$data))){
        #Create something to store the returned data

            datlist[[1]] <- list(data$data)
        }else{
            datlist[[1]] <- data$data
        }
            #cat(str(data))
        # Check to see if there's pagination required
        if ("meta" %in% names(data)){
            #Yes, paginate
            pstart <- as.numeric(data$meta$pagination$page)
            pend <- as.numeric(data$meta$pagination$pages)

            for (p in seq(pstart+1,pend)){    # We've already got the first one

                full_qopts$page <- p
                if(!is.null(client@authtok)){
                    add_headers(.headers = c(Authorization = paste("Bearer", client@authtok, sep=" ")))
                }
                curd <- content(GET(fullurl, config(verbose=Debug), query=full_qopts ))
                datlist[[p]] <- curd$data
                #Check to see if we've pulled enough entries
                if(maxhits > 0){
                    curlen <- sum(sapply(datlist, length))
                    if (curlen > maxhits){
                        break
                    }
                }
            }
        }
        #if(length(datlist) > 1){
        final_data <- unlist(datlist, recursive=F)

        if (usecache && !file.exists(cache_full_fname)){
            #Make sure the directory is created...
        dir.create(dirname(cache_full_fname), recursive = T, showWarnings = client@warnings)
            saveRDS(final_data, file = cache_full_fname)
        }
    }
    final_data
}
