# NOTES:

#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


library(urltools)
library(phyloseq)
library(plyr)
library(reshape2)
library(dplyr)
library(ape)



##Example accessions:
# Analysis assembly: MGYA00379728
# Analysis metagenomic: MGYA00377678
# Analysis amplicon: MGYA00250889



baseurl='https://www.ebi.ac.uk/metagenomics/api/v1'

##Filters possible - this comes from the django source code - would be nice if we could
# look it up.
# These DON'T seem to include all possible attributes ....
# And only some
sample_filters = c('accession','experiment_type','biome_name','lineage','geo_loc_name','latitude_gte','latitude_lte',
                   'longitude_gte','longitude_lte','species','instrument_model','instrument_platform','metadata_key',
                   'metadata_value_gte','metadata_value_lte','metadata_value','environment_material','environment_feature',
                   'study_accession','include')
biome_filters = c('depth_gte','depth_lte')
study_filters = c('accession','biome_name','lineage','centre_name','include')
run_filters=c('accession','experiment_type','biome_name','lineage','species','instrument_platform','instrument_model',
              # 'metadata_key','metadata_value_gte','metadata_value_lte','metadata_value','sample_accession','study_accession',
              'include')
analysis_filters = c('biome_name', 'lineage', 'experiment_type', 'species', 'sample_accession', 'pipeline_version')


#Combined together into a single queriably list
query_filters=list(
  biomes=biome_filters,
  samples=sample_filters,
  studies  = study_filters,
  runs=run_filters
)


### Result table caching

mgnify_memory_cache=list()


#helper function for getting relative paths in the API
#Not everything is implemented here - just what we
#need to get to the download or run areas
#Given an accession x, we want to get the link to get the url for the
#coresponding typeY
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
#' @return \code{char} complete url to access the result. Note this query is not run from here - just the URL is returned
#' @examples
#' cl <- new("mgnify_client")
#' mgnify_get_x_for_y(cl, "MGYS00005126", "studies", "samples")
## @export
mgnify_get_x_for_y <- function(client, x, typeX, typeY, usecache=F){
  #This one's easy - just rearrange the URLs
  #if(typeX=="samples" & typeY %in% c("runs","studies")){
  #  paste( typeX,x,typeY, sep="/")
  #}else if(typeX=="runs" & typeY == "analyses"){
  #  paste( typeX,x,typeY, sep="/")
  #}
  #else{
    #Do it the hard way with a callout
    json_dat = mgnify_retrieve_json(client, paste(typeX, x, sep="/"), usecache = usecache)
    #cat(str(json_dat))
    #tgt_access = json_dat[[1]]$relationships[[typeY]]$data$id
    #tgt_type = json_dat[[1]]$relationships[[typeY]]$data$type
    #paste(tgt_type,tgt_access,sep="/")
    json_dat[[1]]$relationships[[typeY]]$links$related
    #substr(tgt_url, nchar(client@url) + 1, nchar(tgt_url))
  #}
}




#Not exporting this - if people want to they can use the
# rjsonapi functionality. Internally, it takes the "attributes" list
#and converts it into a single row data.frame. For some entries, there is a sublist
#of key/value pairs. metadata_key allows these to be included as columns in the result.
mgnify_attr_list_to_df_row <- function (json, metadata_key=NULL ){

  attrlist=names(json$attributes)

  if (!is.null(metadata_key)){
    baseattrlist=attrlist[!(attrlist %in% c(metadata_key))]
    metaattrlist=json$attributes[[metadata_key]]
    metlist=sapply(metaattrlist, function(x) x$value)
    names(metlist)=sapply(metaattrlist, function(x) x$key)
    df = as.data.frame(t(unlist(c(json$attributes[baseattrlist], metlist))), stringsAsFactors = F)
  }else{
    df = as.data.frame(t(unlist(json["attributes"])), stringsAsFactors = F)
  }
  df$accession <- json$id
  df$acc_type <- json$type

  rownames(df)=df$accession
  df
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
  study_df <- mgnify_attr_list_to_df_row(study_met[[1]])

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
    error=function(x) warning("Error finding biome entry")
  )

  full_df
}



#Internal functions to parse the attributes/hierarchy list into a data.frame
mgnify_parse_tax <- function(json){
  df <- as.data.frame(c(json$attributes["count"], unlist(json$attributes$hierarchy)), stringsAsFactors = F)
  df$index_id <- json$attributes$lineage
  df

}
mgnify_parse_func <- function(json){
  df <- as.data.frame(json$attributes, stringsAsFactors = F)
  df$index_id <- json$attributes$accession
  df
}

#Which parser do you use for which type of output?
analyses_results_type_parsers <- list(taxonomy=mgnify_parse_tax,`taxonomy-itsonedb` = mgnify_parse_tax, `go-slim`=mgnify_parse_func,
                                      `taxonomy-itsunite`=mgnify_parse_tax, `taxonomy-ssu`=mgnify_parse_tax,
                                      `taxonomy-lsu`=mgnify_parse_tax,`antismash-gene-clusters`=mgnify_parse_func,
                                      `go-terms`=mgnify_parse_func, `interpro-identifiers`=mgnify_parse_func)


#this maps the json attribute name for retrievelist to the "description -> label" attribute in the study downloads section
analyses_results_bulk_file_names <- list(
                                    #taxonomy="Taxonomic assignments SSU",
                                    `taxonomy-itsonedb` = "Taxonomic assignments ITS",
                                    `go-slim`="GO slim annotation",
                                    `taxonomy-itsunite`="Taxonomic assignments Unite", `taxonomy-ssu`="Taxonomic assignments SSU",
                                    `taxonomy-lsu`="Taxonomic assignments LSU",`antismash-gene-clusters`=mgnify_parse_func,
                                    `go-terms`="Complete GO annotation", `interpro-identifiers`="InterPro matches",
                                    `phylo-tax-ssu`="Phylum level taxonomies SSU",`phylo-tax-lsu`="Phylum level taxonomies LSU" )


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
          urltools::parameters(data_url) <- NULL

          #build the cache filename
          fname=tail(strsplit(data_url, '/')[[1]], n=1)

          #At this point we might have alread got the data we want loaded. Check the memory cache object

          if(client@use_memcache & (cur_type %in% names(mgnify_memory_cache)) & (mgnify_memory_cache[cur_type]["fname"] == fname)){
            tmp_df <- mgnify_memory_cache[cur_type][["data"]]
          }else{
            #Nope - gonna have to load it up from disk or grab it from t'interweb
            data_path = paste(downloadDIR, fname, sep="/")
            if (! file.exists(data_path)){#} | !use_downloads ){
              httr::GET(data_url, httr::write_disk(data_path, overwrite = T ))
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
          i=1
          #tmp_df <- read.csv('~/.MGnify_cache/tsv/ERP108138_IPR_abundances_v4.1.tsv', sep="\t", header = T, stringsAsFactors = F)
          while(any(is.na(suppressWarnings(as.numeric(tmp_df[,i] ))))){
            i=i+1
          }
          i=i-1

          #also need the column name for this particular analysis... As far as I can see they could be either assembly IDs or run ids. FFS.
          #Assuming that both assembly and run won't be present...:
          if("assembly_accession" %in% colnames(metadata_df)){
            accession=metadata_df$assembly_accession[[1]]
          }else if("run_accession" %in% colnames(metadata_df)){
            accession=metadata_df$run_accession[[1]]
          }
        #cat(accession)
        #  cat(str(head(tmp_df[1:5,1:5])))
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
    names(parsed_results) <-  names(analyses_results_bulk_file_names)[match(unlist(lapply(available_downloads_json,
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
    parsed_results = sapply(names(all_results), function(x){
      all_json <- all_results[[x]]
      if(! is.null(all_json)){
        res_df <- do.call(dplyr::bind_rows, lapply(all_json,analyses_results_type_parsers[[x]] ))
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


#UPDATE ME SO THAT TREES (if available) GET GRABBED AS WELL!!!
# Not exported - get a single biom file and convert it to a phyloseq object.
mgnify_get_single_analysis_phyloseq <- function(client=NULL, accession, usecache=T, downloadDIR=NULL, tax_SU="SSU", get_tree=FALSE){
  metadata_df <- mgnify_get_single_analysis_metadata(client, accession, usecache=usecache)

  analysis_data <-  mgnify_retrieve_json(client, paste("analyses",accession,sep="/"), usecache = usecache)
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
    downloadDIR=paste(client@cache_dir,"biom_files",sep="/")
    dir.create(downloadDIR, recursive = T, showWarnings = client@warnings)
  }
  #Clear out any ?params after the main location - don't need them for this
  urltools::parameters(biom_url) <- NULL

  fname=tail(strsplit(biom_url, '/')[[1]], n=1)
  biom_path = paste(downloadDIR, fname, sep="/")
  if (! file.exists(biom_path)){#} | !use_downloads ){
    httr::GET(biom_url, httr::write_disk(biom_path, overwrite = T))
  }
  #Load in the phlyloseq object
  psobj <- phyloseq::import_biom(biom_path)
  #Need to check if the taxonomy was parsed correctly - depending on the pipeline it may need a bit of help:
  if (ncol(phyloseq::tax_table(psobj)) == 1){
    psobj <- phyloseq::import_biom(biom_path, parseFunction = phyloseq::parse_taxonomy_qiime)
  }
  if(! "Kingdom" %in% names(phyloseq::tax_table(psobj))){
    psobj <- phyloseq::import_biom(biom_path, parseFunction = phyloseq::parse_taxonomy_greengenes)
  }

  #The biom files have a single column of unknown name - maybe it's the original sample name?.
  # It's rewritten as sample_run_analysis accession, with
  # the original value stored in the sample_data (just in case it's needed later)
  orig_samp_name <- phyloseq::sample_names(psobj)[[1]]
  newsampname <- rownames(metadata_df)[1]
  metadata_df[1,"orig_samp_name"] <- orig_samp_name
  phyloseq::sample_names(psobj) <- newsampname
  phyloseq::sample_data(psobj) <- metadata_df

  #Finally, do we want to the phylogenetic tree? If so, is it there?
  if(get_tree){
    #is there a tree?
    tvec = grepl('Phylogenetic tree', sapply(analysis_downloads, function(x) x$attributes$`description`$label))
    if(any(tvec)){
      tree_url = analysis_downloads[tvec][[1]]$links$self
      #Clear out any ?params after the main location - don't need them for this
      urltools::parameters(tree_url) <- NULL

      fname=tail(strsplit(tree_url, '/')[[1]], n=1)
      tree_path = paste(downloadDIR, fname, sep="/")
      if (! file.exists(tree_path)){#} | !use_downloads ){
        httr::GET(tree_url, httr::write_disk(tree_path, overwrite = T ))
      }
    }
    phylo_tree = ape::read.tree(tree_path)
    phyloseq::phy_tree(psobj) <- phylo_tree
  }
  psobj
}



##' @exportClass mgnify_client
mgnify_client <- setClass("mgnify_client",
                          slots=list(url = "character", authtok = "character",
                                     cache_dir="character", warnings="logical",
                                     use_memcache="logical", memcache="list"),
                          prototype = list(url=baseurl, authtok=NULL, cache_dir=NULL, use_memcache=FALSE, memcache=list()))

#Contructor to allow logging in with username/password
#' Instantiate the MGnifyR client object
#'
#' All functions in the MGnifyR package take a \code{mgnify_client} object as their first argument. While not essential
#' to querying the raw MGnify API (which is exposed as relative standard JSONAPI), the object allows the simple handling of both
#' user authentication and access to private data, and local on-disk caching of results.
#'
#' @param username optional username to authenticate.
#' @param password optional password for authentication.
#' @param usecache whether to enable on-disk caching of results during this session. In most use cases should be TRUE.
#' @param cache_dir specifies a folder to contain the local cache. If NULL, and usecache is TRUE, the new subdirectory \code{.MGnifyR_cache}
#'  in the current working directory will be used. Note that cached files are persistent, so the cache directory may be reused between sessions,
#'  taking advantage of previously downloaded results. The directory will be created if it doesn't exist already.
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
    url=baseurl
  }

  authtok=NA_character_

  #Check to see if we're goint to try and get an authentication token:
  if (!is.null(username) & !is.null(password)){
    r = httr::POST(paste(url, "utils/token/obtain", sep="/"),
                   body=list(username=username, password=password),
                   encode="json")
    cont = httr::content(r)
    if ("data" %in% names(cont)){
      authtok = cont$data$token
    }
    else{
      stop("Failed to authenticate")
    }
  }
  #Assume we're not using it
  cachepath=NA_character_
  if(usecache){
    if (is.null(cache_dir) ){
      cachepath=paste(getwd(),'.MGnifyR_cache',sep="/")
    }else{
      cachepath=cache_dir
    }
    #Make it if needed - assume the user is sensible and the path will work...
    dir.create(cachepath,showWarnings = F)
  }

  #Return the final object
  new("mgnify_client", url=url, authtok=authtok, cache_dir = cachepath, warnings=warnings, memcache=list(), use_memcache=use_memcache)
}


#Internal function to actually perform the http request. Build up the URL then issues
#a GET, parsing the returned JSON into a nested list (uses \code{jsonlite} internally?)
#Previously cached results may be retrieved from disk without resorting to calling the MGnify server.

#'Low level MGnify API handler
#'
#'\code{mgnify_retrieve_json} deals with handles the actual HTTP GET calls for the MGnifyR package, handling API pagination,
#'local result caching, and  authentication cookies for access
#'to restricted or pre-release datasets.Although principally intended for internal MGnifyR use , it's exported for direct invocation.
#'Generally though it's not recommended for use by users.
#'
#'@param client MGnifyR client
#'@param path top level search point for the query. One of \code{biomes}, \code{samples}, \code{runs} etc. Basically includes
#'all parts of the URL between the base API url and the parameter specifications
#'@param complete_url \emph{complete} url to search, usuaally retrieved from previous query in the "related" section.
#'@param qopts named list or vector containing options/filters to be URL encoded and appended to query as key/value pairs
#'@param maxhits Maxmium number of data entries to return. The actual number of hits returned may be higher than this value,
#'as this parameter only clamps after each full page is processed. Set to <=0 to disable - i.e. retrieve all items.
#'@param usecache Should successful queries be cached on disk locally? There are unresolved questions about whether this is
#'a sensible thing to do, but it remains as an option. It probably makes sense for single accession grabs, but not for
#'(filtered) queries - which are liable to change as new data is added to MGnify. Also caching only works for the first page.
#'@param Debug Should we print out lots of information while doing the grabbing?
#'@return \code{list} of results after pagination is dealt with.
#'@export
  mgnify_retrieve_json <- function(client, path="biomes", complete_url=NULL, qopts=NULL,
                                   maxhits=200, usecache = F, Debug=F){
  #Set up the base url
  #Are we using internal paths?
  if (is.null(complete_url)){
    fullurl = paste(client@url, path, sep="/")
  }
  #Or direct links from e.g. a "related" section
  else{
    #Set the full url, but clean off any existing parameters (page, format etc) as they'll be added back later:
    fullurl = complete_url
    urltools::parameters(fullurl) <- NULL
    path = substr(fullurl, nchar(client@url) + 2, nchar(fullurl))
  }

  #cat(fullurl)

  #convert to csv if filters are lists.
  #This doesn't check if they ~can~ be searched for in the API,
  #which is an issue since no error is returned by the JSON if the search
  #is invalid - we only get a result as if no query was present...
  tmpqopts = lapply(qopts,function(x) paste(x,collapse = ','))

  #Include the json and page position options
  #full_qopts = as.list(c(format="json", tmpqopts, page=1))
  full_qopts = as.list(c(format="json", tmpqopts))
  #Build up the cache name anyway - even if it's not ultimately used:
  fname_list = c(path, names(unlist(full_qopts)), unlist(full_qopts))
  cache_fname = paste(fname_list,collapse = "_")
  cache_full_fname = paste(client@cache_dir, '/', cache_fname, '.RDS', sep="")

  # Do we want to try and use a cache to speed things up?
  if(usecache & file.exists(cache_full_fname)){
      final_data = readRDS(cache_full_fname)
  }else{

    #Authorization: Bearer <your_token>
    if(!is.null(client@authtok)){
      httr::add_headers(.headers = c(Authorization = paste("Bearer", client@authtok, sep=" ")))
    }
    res = httr::GET(url=fullurl, httr::config(verbose=Debug), query=full_qopts )
    data <-httr::content(res)

    #At this point, data$data is either a list of lists or a single named list. If it's a single entry, it needs embedding in
    #a list for consistency downstream
    #datlist is built up as a list of pages, where each entry must be another list. Thus, on the first page,
    #
    datlist=list()
    if (!is.null(names(data$data))){
    #Create something to store the returned data

      datlist[[1]] = list(data$data)
    }else{
      datlist[[1]] = data$data
    }
      #cat(str(data))
    # Check to see if there's pagination required
    if ("meta" %in% names(data)){
      #Yes, paginate
      pstart = as.numeric(data$meta$pagination$page)
      pend   = as.numeric(data$meta$pagination$pages)

      for (p in seq(pstart+1,pend)){  # We've already got the first one

        full_qopts$page=p
        if(!is.null(client@authtok)){
          httr::add_headers(.headers = c(Authorization = paste("Bearer", client@authtok, sep=" ")))
        }
        curd = httr::content(httr::GET(fullurl, httr::config(verbose=Debug), query=full_qopts ))
        datlist[[p]] = curd$data
        #Check to see if we've pulled enough entries
        if(maxhits > 0){
          curlen=sum(sapply(datlist, length))
          if (curlen > maxhits){
            break
          }
        }
      }
    }
    #if(length(datlist) > 1){
    final_data <- unlist(datlist, recursive=F)
    #}else{
    #  final_data <- datlist
    #}
    #}
    #else{
    #  final_data <- datlist
    if (usecache && !file.exists(cache_full_fname)){
      #Make sure the directory is created...
    dir.create(dirname(cache_full_fname), recursive = T, showWarnings = client@warnings)
      saveRDS(final_data, file = cache_full_fname)
    }
  }
  final_data
}


#' Listing files available for download
#'
#' \code{mgnify_get_downloads} is a wrapper function allowing easy enumeration of downloads available for a given
#' accession (or list thereof). Returns a single data.frame containing all available downloads and associated metadata,
#' including the url location and description. This can then be filtered to extract the urls of interest, before actually
#' retrieving the files using \code{mgnify_download}
#'
#'@param client valid MGnify client object
#'@param accessions list of accessions to query
#'@param accession_type one of \code{analysis},\code{samples},\code{studies},\code{assembly},\code{genome} or \code{run}
#'@param usecache whether to use the on-disk cache to speed up queries (default T)
#'@return \code{data.frame} containing all discovered downloads. If multiple \code{accessions} are queried, the \code{accessions} column
#' may to filter the results - since rownames are not set (and wouldn;'t make sense as each query will return multiple items)

#'@examples
#' #Make a client ibject
#' mg <- mgnify_client(cache_dir="/tmp/mgcache")
#' #create a vector of accession ids - these happen to be \code{analysis} accessions
#' accession_vect <- c("MGYA00563876", "MGYA00563877", "MGYA00563878", "MGYA00563879", "MGYA00563880" )
#' downloads <- mgnify_get_downloads(mg, accession_vect, "analyses")
#'@export

mgnify_get_download_urls <- function(client, accessions, accession_type, usecache=T){
  results <- plyr::llply(accessions, function(x){
    download_list  <- mgnify_retrieve_json(client, paste(accession_type,x,"downloads", sep="/"), usecache = usecache)
    df <- do.call(rbind.fill,lapply(download_list, function(x) as.data.frame(x,stringsAsFactors=F)))
    df$accession <- x
    df$accession_type <- accession_type
    #for convenience, rename the "self" column to "download_url" - which is what it actually is...
    colnames(df)[colnames(df) == 'self'] <- 'download_url'
    #finally, strip off any options from the url - they sometimes seem to get format=json stuck on the end
    urls <- df$download_url
    urltools::parameters(urls) <- NULL
    df$download_url <- urls
    df
  }, .progress="text")
  do.call(rbind.fill, results)
}


#' Download arbitray files from MGnify, including processed reads and identified protein sequences.
#'
#' \code{mgnify_download} is a convenient wrapper round generic the url downloading functionality in R, taking care of things like local
#' caching and authentication. By default, \code{mgnify_download}
#' @param client MGnify client object
#' @param url The url of the file we wish to download
#' @param target_filename An optional local filename to use for saving the file. If NULL (default), MGNify local cache settings will be used.
#' If the file is intended to be processed in a seperate program, it may be sensible to provide a meaningful \code{target_filename}, rather than having to hunt
#' through the cache folders. If \code{target_filename} is NULL \emph{and} \code{usecache} is \code{FALSE}, the \code{read_func} parameter must be supplied or the file
#' will be downloaded and then deleted.
#' @param read_func An optional function name to process the downloaded file and return the results, rather than relying on post processing. The primary use=case for
#'this parameter is when local disk space is limited and downloaded files can be quickly processed and discarded. The function should take a single parameter,
#'the downloaded filename, and may return any valid R object.
#' @param usecache whether to enable the default MGnifyR caching mechanism. File locations are overridden if \code{target_filename} is supplied.
#' @param Debug whether to enable debug output of the HTTP call - only useful for development.
#' @return Either the local filename of the downloaded file, be it either the location in the MGNifyR cache, or target_filename. If \code{read_func} is used, its result
#' will be returned.
#' @examples
#' #Make a client object
#' mg <- mgnify_client(cache_dir="/tmp/mgcache")
#' #create a vector of accession ids - these happen to be \code{analysis} accessions
#' accession_vect <- c("MGYA00563876", "MGYA00563877", "MGYA00563878", "MGYA00563879", "MGYA00563880" )
#' downloads <- mgnify_get_downloads(mg, accession_vect, "analyses")
#'
#' #Filter to find the urls of 16S encoding sequences
#' url_list <- downloads[downloads$attributes.description.label == "Contigs encoding SSU rRNA","download_url"]
#'
#' #Example 1:
#' #Download the first file
#' supplied_filename = mgnify_download(mg, url_list[[1]], target_filename="SSU_file.fasta.gz")
#'
#'
#' #Example 2:
#' #Just use local caching
#' cached_filename = mgnify_download(mg, url_list[[2]])
#'
#' #Example 3:
#' #Using read_func to open the reads with readDNAStringSet from \code{biostrings}. Without retaining on disk
#' dna_seqs <- mgnify_download(mg, url_list[[3]], read_func=readDNAStringSet, usecache=F)
#'
#' @export
mgnify_download <- function(client, url, target_filename=NULL, read_func=NULL, usecache=TRUE, Debug=FALSE){
  #Set up filenames for storing the data
  ftgt=NULL
  if (! is.null(target_filename)){
    file_tgt = target_filename
  }else if(usecache == TRUE){
    #Build a filename out of the url, including the full paths. Annoying, but some downloads (e.g. genome results) are just names like
    # core_genes.fa , which would break the caching.
    cachetgt = gsub(paste(client@url,'/',sep=""), '', url)
    #Make sure the direcory exists

    cache_full_name = paste(client@cache_dir, cachetgt, sep="/")
    dir.create(dirname(cache_full_name), recursive = T, showWarnings = client@warnings)


    file_tgt = cache_full_name
  } else{
    file_tgt = tempfile()[[1]]
  }

  #Only get the data if it's not already on disk
  if(!(usecache & file.exists(file_tgt))){

    if(!is.null(client@authtok)){
      httr::add_headers(.headers = c(Authorization = paste("Bearer", client@authtok, sep=" ")))
    }
    #If there's an error we need to make sure the cache file isn't written - by default it seems it is.
    tryCatch(
    curd = httr::content(httr::GET(url, httr::write_disk(file_tgt, overwrite = T)))
    , error=function(x){
      unlink(file_tgt)
      print(paste("Error retrieving file",file_tgt))
      print(paste("Error:",x))
      stop()
    })
  }

  if (is.null(read_func)){
    result = file_tgt
  } else{
    result = read_func(file_tgt)
  }

  if (is.null(target_filename) & !usecache){
    #Need to clear out the temporary file
    unlink(file_tgt)
  }
  result
}




#'#' Search MGnify database for studies, samples and runs
#'
#' \code{mgnify_query} is a flexible query function, harnessing the "full" power of the JSONAPI MGnify
#' search filters. Search results may be filtered by metadata value, associated study/sample/analyese etc. Details
#' of the capabilities may be found [here](https://emg-docs.readthedocs.io/en/latest/api.html#customising-queries). Currently,
#' the following filters are available (based on examination of the Python source code):
#'   \itemize{
#'     \item{\strong{Studies} : accession, biome_name, lineage, centre_name}
#'      \item{\strong{Samples} : accession, experiment_type, biome_name,
#'   lineage, geo_loc_name, latitude_gte, latitude_lte,
#'   longitude_gte, longitude_lte, species, instrument_model, instrument_platform,
#'    metadata_key, metadata_value_gte, metadata_value_lte, metadata_value,
#'    environment_material, environment_feature, study_accession}
#'    \item{\strong{Runs} accession, experiment_type, biome_name, lineage, species,
#'      instrument_platform, instrument_model}
#'    }
#'    Unfortunately (from testing) it appears that some of these filters don't work as expected, so it is important
#'    to check the results returned match up with what's expected. Even more unfortunately if there's an error in the
#'    parameter specification, the query will run as if no filter parameters were present at all. Thus the
#'    result will appear superficially correct but will infact correspond to something completely different. This beahviour
#'    will hopefully be fixed in future incarnations of the API, but for now users should double check returned
#'    values.
#'
#'    It is currently not possible to combine queries of the same type in a single call (for example to search for samples
#'    \emph{between} latitude). However, it is possible to run multiple queries and combine the results using set operations in R to get the
#'    desired behaviour.
#'
#'
#'
#' @param mgnify_client Client instance
#' @param qtype Type of objects to query. One of \code{studies},\code{samples},\code{runs} or
#' \code{analyses}
#' @param accession Either a single known MGnify accession identifier (of type \code{qtype}), or a list/vector
#' of accessions to query. Note that multiple values only work for samples, runs and assemblies ... not sure why.
#' @param asDataFrame Boolean flag to choose whether to return the results as a data.frame or leave as a nested list. In
#' most cases, \code{asDataFrame = TRUE} will make the most sense.
#' @param maxhits determines the maximum number of results to return. The actual number of results will actually be higher than \code{maxhits},
#' as clipping only occurs on pagination page boundaries. To disable the limit, set \code{maxhits} < 0
#' @param usecache Whether to cache the result - and reuse any existing cache entry instead of issuing a
#' new callout. In generl the use of caching for queries is discouraged, as new data is being uploaded to MGnify
#' all the time, which might potentially be missed. However, for some purposes (such as analysis reproducibility)
#' caching makes sense.
#' @param ... Remaining parameter key/value pairs may be supplied to filter the returned values. Available options differ
#' between \code{qtypes}.See discussion above for details.
#' @return  A nested list or data.frame (depending on \code{asDataFrame}) containing the results of the query.
#' @examples
#' mg <- mgnify_client(cache_dir="/tmp/mgcache")
#'
#' #Get a list of studies from the Agricultural Wastewater :
#' agwaste_studies <- mgnify_query(mg, "studies", biome_name="Agricultural wastewater")
#'
#' #Get all samples from a particular study
#' samps <- mgnify_query(mg, "samples", study_accession="MGYS00004521")
#'
#' #Search for all polar samples
#' samps_np <- mgnify_query(mg, "samples", latitude_gte=66, maxhits=-1)
#' samps_sp <- mgnify_query(mg, "samples", latitude_lte=-66, maxhits=-1)
#' samps_polar <- rbind(samps_np, samps_sp)
#'
#' @export
mgnify_query <- function(client, qtype="samples", accession=NULL, asDataFrame=T, maxhits=200, usecache=F, ...){
  #Need to get around the lazy expansion in R in order to get a list
  a=accession
  arglist = as.list(match.call())[-1] # drop off the first entry, which is the name of the function

  arglist$accession=a

  #Filter the query options such that
  qopt_list = arglist[names(arglist) %in% query_filters[[qtype]]]
  non_qopts = arglist[!(names(arglist) %in% c(c("asDataFrame","qtype","client", "maxhits"),query_filters[[qtype]]))]

  all_query_params = unlist(list(c(list(client=client, maxhits=maxhits, path=qtype, usecache=usecache, qopts=qopt_list))), recursive = F)

  result = do.call("mgnify_retrieve_json", all_query_params)

  #Rename entries by accession
  id_list = lapply(result, function(x) x$id)
  names(result) = id_list

  if(asDataFrame){
    #Because metadata might not match across studies, the full dataframe is built by first building per-sample dataframes,
    # then using rbind.fill from plyr to combine. For ~most~ use cases the number of empty columns will hopefully
    # be minimal... because who's going to want cross study grabbing (?)
    dflist = lapply(result, function(r){
      df2 <- mgnify_attr_list_to_df_row(json = r, metadata_key = "sample-metadata")

      #Currently this is a bit hacky - assumes the study only has one biome, and sample only one study etc.
      for(rn in names(r$relationships)){
        tryCatch({
          df2[rn] = as.list(r$relationships[[rn]]$data)[[1]]$id
        },
        error=function(x)NULL)
      }
      df2$type = r$type
      rownames(df2)=df2$accession
      df2

    }
    )
    tryCatch(
      dplyr::bind_rows(dflist),
      error=function(e) dflist
    )
  }else{
    result
  }

}


#' Look up analysis accession IDs for one or more study accessions
#'
#' \code{mgnify_analyses_from_studies} Retrieve Analysis accession IDs associated with the supplied Study accession
#'
#' Helper function to get all analyses associated with the given studies.
#'
#' @param client \code{mgnify_client} instance
#' @param accession Single study accession id, or vector/list of accessions for which to retrieve Analyses ids
#' @param usecache Flag to determine whether to re-use/store data on disk, rather than query the server.
#' @return vector of Analysis accession ids
#' @examples
#' #Retrieve all analysis ids from studies MGYS00005058, MGYS00005058 and MGYS00005058
#' result <- mgnify_analyses_from_studies(myclient, c("MGYS00005058", "MGYS00005058" and "MGYS00005058"))
#'
#' @export
mgnify_analyses_from_studies <- function(client, accession, usecache=T){
  analyses_accessions <- plyr::llply(as.list(accession), function(x){
    accurl <- mgnify_get_x_for_y(client, x, "studies","analyses", usecache = usecache )
    jsondat <- mgnify_retrieve_json(client, complete_url = accurl, usecache = usecache, maxhits = -1)
    #Just need the accession ID
    lapply(jsondat, function(x) x$id)
  }, .progress="text")
  unlist(analyses_accessions)
}


#' Look up analysis accession IDs for one or more sample accessions
#'
#' \code{mgnify_analyses_from_samples} Retrieve Analysis accession IDs associated with the supplied Sample accession
#'
#' Helper function to get all analyses associated with the given samples.
#'
#' @param client \code{mgnify_client} instance
#' @param accession Single sample accession id, or vector/list of accessions for which to retrieve Analyses ids
#' @param usecache Flag to determine whether to re-use/store data on disk, rather than query the server.
#' @return vector of associated Analysis accession ids
#' @examples
#' #Retrieve all analysis ids from samples
#' result <- mgnify_analyses_from_samples(myclient, c("SRS4392730", "SRS4392743"))
#'
#' @export
mgnify_analyses_from_samples <- function(client, accession, usecache=T){
  #analyses_accessions <- sapply(as.list(accession), function(x){
  analyses_accessions <- plyr::llply(as.list(accession), function(x){
    accurl <- mgnify_get_x_for_y(client, x, "samples","analyses", usecache = usecache )
    #For some reason, it appears you "sometimes" have to go from study to runs to analyses. Need
    #to query this with the API people...
    if(is.null(accurl)){
      runurl <- mgnify_get_x_for_y(client, x, "samples","runs", usecache = usecache )
      jsondat <- mgnify_retrieve_json(client, complete_url = runurl, usecache = usecache)
      run_accs <- lapply(jsondat, function(y) y$id)
      a_access <- sapply(as.list(run_accs), function(z){
        accurl <- mgnify_get_x_for_y(client, z, "runs","analyses", usecache = usecache )
        jsondat <- mgnify_retrieve_json(client, complete_url = accurl, usecache = usecache)
        lapply(jsondat, function(x) x$id)
      })
      unlist(a_access)
    }else{
      jsondat <- mgnify_retrieve_json(client, complete_url = accurl, usecache = usecache)
      #Just need the accession ID
      lapply(jsondat, function(x) x$id)
    }}, .progress="text")
  unlist(analyses_accessions)
}




#' Get all Study, Sample and Analysis metadata for the supplied analyses accessions
#'
#' \code{mgnify_get_analyses_metadata} retrieves all associated Study, Sample and Analysis metadata attributes
#' a list of Analyses accessions (determined from \code{mgnify_analyses_from_x})
#'
#' @param client \code{mgnify_client} instance
#' @param accessions Single value or list/vector of Anlysis accessions to retrieve data for
#' @param usecache Whether to use the disk based cache.
#' @return \code{data.frame} of metadta for each analysis in the \code{accession} list.
#' @examples
#'
#' @export
mgnify_get_analyses_metadata <- function(client, accessions, usecache=T){
  reslist <- plyr::llply(as.list(accessions), function(x) mgnify_get_single_analysis_metadata(client, x, usecache = usecache),
                            .progress = "text")
  df <- do.call(dplyr::bind_rows,reslist)
  rownames(df) <- accessions
  df
}




#' Retrieve OTU tables for all specified accessions and build a\code{phyloseq} object
#'
#' \code{mgnify_get_analyses_phyloseq} retrieves all associated Study, Sample and Analysis metadata attributes,
#' along with all OTU tables (of a given taxonomic type), and merges them together to build a \code{phyloseq}
#' object with \code{otu_table}, \code{tax_table} and \code{sample_data} objects.
#' #'
#' @param client \code{mgnify_client} instance
#' @param accessions Single value or list/vector of Analysis accessions to retrieve data for.
#' @param usecache Whether to use the disk based cache.
#' @param returnLists Flag to determine whether to merge a per-analysis \code{phyloseq} objects into one, or
#' return a list of single-sample objects. Since perculiarities in OTU format may lead to failed merges, returning
#' a list of single objects may aid debugging. In most cases though the desired behaviour is to return a single
#' object (returnLists = F)
#' @param tax_SU Which taxa subunit results should be selected? Currently, taxonomy assignments in the
#' MGnify pipelines rely on rRNA matches to existing databases (GreenGenes and SILVA), with later
#' pipelines checking both the SSU and LSU portions of the rRNA sequence. \code{tax_SU} allows the
#' selection of either the Small subunit (SSU) or Large subunit results in the final \code{phyloseq} object.
#' Older pipeline versions do not report results for both subunits,
#' and thus for some accessions this value will have no effect.
#' @param get_tree Flag to control whether to include available phylogenetic trees in the phyloseq object. At present this option is of limited use as
#' the trees available are specific to each accession, holding only a subset of the leaves and branches of the full canonical tree used for the pipeline. This means that
#' \code{phyloseq} is unable to merge the trees, and therefore fails to build a final combined object. Setting \code{returnLists} to TRUE allows the results to be returned successfully,
#' albeit not in a single object.
#' @return
#' @examples
#'
#'
#' @export
mgnify_get_analyses_phyloseq <- function(client = NULL, accessions, usecache=T,
                                         returnLists=F, tax_SU = "SSU",
                                         get_tree=FALSE){
  #Some biom files don't import - so need a try/catch
  ps_list <- plyr::llply(accessions, function(x) {
    tryCatch(
        mgnify_get_single_analysis_phyloseq(client, x, usecache = usecache, tax_SU = tax_SU, get_tree = get_tree), error=function(x) NULL)
    }, .progress = "text")

  #The sample_data has been corrupted by doing the merge (names get messed up and duplicated), so just regrab it with another lapply/rbind
  samp_dat <- lapply(accessions, function(x) mgnify_get_single_analysis_metadata(client, x, usecache = usecache ))
  if (returnLists){
    list(phyloseq_objects=ps_list, sample_metadata = samp_dat)
  }else{

    #first of all, check to see that if we wanted them, we got trees for ALL the phyloseq objects.
    #If trees are present in any of the phyloseq objects during merging, then any OTUs not in a tree
    #(e.g. if any phyloseq objects do NOT contain a tree) will not be included in the merged output.

    if(get_tree){
      if (any(is.na(lapply(ps_list, function(x) x@phy_tree)))){
      warning("Phylogenetic tree retrieval was requested but some of the analyses do not include phylogenetic trees. Results should be used with caution.")
      }
    }

    #This is too slow for large datasets
    #full_ps <- do.call(phyloseq::merge_phyloseq, ps_list)
    #so:
    #a divide-and-conquer approach to merge_phyloseq seems to work best, hence the following
    #code which splits the full list into sublists and merges them seperately, then repeats until all are joined.
    curlist=ps_list
    while(length(curlist) > 1){
      #Lists of length 10 seem to work well
      sublist=split(curlist, seq_along(curlist) %/% 10)
      curlist <- lapply(sublist, function(x){
        do.call(phyloseq::merge_phyloseq,x)
      })
    }
    #By this point curlist isn't a list, it's a phyloseq object...
    full_ps <- curlist[[1]]

    sample_metadata_df <- do.call(dplyr::bind_rows, samp_dat)
    rownames(sample_metadata_df) <- sample_metadata_df$analysis_accession
    phyloseq::sample_data(full_ps) <- sample_metadata_df
    full_ps
  }
}


#' Get functional and taxonomic information for a list of accessions
#'
#' Given a set of analysis accessions and collection of annotation types, \code{mgnify_get_analyses_results} queries the MGNify API
#' and returns the results, by default merging the results into multi-accession data.frames
#'
#' @param client a valid \code{mgnify_client} object
#' @param accessions list or vector of accessions to return results for
#' @param retrievelist list or vector of functional analysis types to retrieve, or "all" to get all available results. The current list of available
#' types can be found using \code{ names(MGnifyR::analyses_results_type_parsers)}. Note that not depending on the particular analysis type, puipeline
#' version etc., not all functional results will be available.
#' @param compact_results optional parameter to return a named list (one entry per element in \code{retrievelist}) of data.frames, with each data.frame
#' containing results for all requested accessions. If \code{FALSE}, \code{mgnify_get_analyses_results} returns a lists of lists, each element consiting of
#' results for a single accession.
#' @param usecache Whether to use the MGnify local caching system to speed up searching. It is highly recommended that this be enabled (default=TRUE)
#' @param bulk_dl should MGnifyR attempt to speed things up by downloading relevant studies TSV results and only extracting the required columns, rather than using
#' the JSONAPI interface. When getting results where multiple accessions share the same study, this option may result in significantly faster processing. However, there
#' appear to be (quite a few) cases in the database where the TSV result columns do NOT match the expected accession names. This will hopefully be fixed in the future, but for
#' now \code{bulk_dl} defaults to FALSE. When it does work, it can be orders of magnitude more efficient.
#' @return
#' @examples
#'
#'@export
mgnify_get_analyses_results <- function(client=NULL, accessions, retrievelist=c(), compact_results=T, usecache = T, bulk_dl = F){
  if(length(retrievelist) == 1 && retrievelist == "all"){
    retrievelist = names(analyses_results_type_parsers)
  }
  results_as_lists <- plyr::llply(accessions,
                                  function(x) mgnify_get_single_analysis_results(
                                    client, x,
                                    usecache = usecache,
                                    retrievelist = retrievelist, bulk_files = bulk_dl),
                                  .progress = "text")
  names(results_as_lists) <- accessions

  if(!compact_results){
    results_as_lists
  }else{
    #Compact the result type dataframes into a single instance. Per accession counts in each column.
    all_results <- plyr::llply(retrievelist, function(y){
      tryCatch({
        r = lapply(results_as_lists, function(x){
          df <- as.data.frame(x[[y]])
          df
        })
        longform <- dplyr::bind_rows(r, .id = "analysis")
        cn <- colnames(longform)
        extras <- cn[!(cn %in% c("count","index_id", "analysis"))]
        final_df <- reshape2::dcast(longform, as.formula(paste(paste(extras,  collapse = " + "), " ~ analysis")), value.var = "count", fun.aggregate = sum)
        final_df}, error=function(x) NULL)
    })
  }
  names(all_results) <- retrievelist
  all_results
}




