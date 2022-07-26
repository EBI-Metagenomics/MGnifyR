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

                        if(usecache & client@clear_cache){
                            message(paste("clear_cache is TRUE: deleting ", data_path, sep=""))
                            tryCatch(unlink(data_path), error=warning)
                        }

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
