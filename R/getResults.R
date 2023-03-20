#' Get functional and/or taxonomic information for a list of accessions
#'
#' @details
#' Given a set of analysis accessions and collection of annotation types,
#' the function queries the MGNify API and returns the results.
#'
#' @param x A \code{MgnifyClient} object.
#'
#' @param accession A single character value or a vector of character values
#' specifying accession IDs to return results for.
#'
#' @param output A single character value specifying the format of an output.
#' Must be one of the following options: "TreeSE", "list", or "phyloseq".
#' (By default: \code{output = "TreeSE"})
#'
#' @param get.taxa A boolean value specifying whether to retrieve metagenomic
#' data. (By default: \code{get.taxa = TRUE})
#'
#' @param get.func A boolean value or a single character value or a vector
#' character values specifying functional analysis types to retrieve. If
#' \code{get.func = TRUE}, all available functional datatypes are retrieved,
#' and if \code{FALSE}, functional data is not retrieved. The current list of
#' available types is "antismash-gene-clusters", "go-slim", "go-terms",
#' "interpro-identifiers", "taxonomy", "taxonomy-itsonedb",
#' "taxonomy-itsunite", "taxonomy-lsu", and "taxonomy-ssu". Note that
#' depending on the particular analysis type, pipeline version etc., not all
#' functional results will be available. (By default: \code{get.func = FALSE})
#'
#' @param use.cache A single boolean value specifying whether to use the
#' MGnify local caching system to speed up searching. It is highly
#' recommended that this be enabled (By default: \code{use.cache = TRUE})
#'
#' @param ... optional arguments:
#' \itemize{
#'   \item{taxa.su}{ A single character value specifying which taxa subunit
#'   results should be selected? Currently, taxonomy assignments in the
#'   MGnify pipelines rely on rRNA matches to existing databases
#'   (GreenGenes and SILVA), with later pipelines checking both the SSU and
#'   LSU portions of the rRNA sequence. \code{taxa.su} allows then selection
#'   of either the Small subunit (SSU) or Large subunit results
#'   in the final \code{TreeSummarizedExperiment} object. Older pipeline
#'   versions do not report results for both subunits, and thus for some
#'   accessions this value will have no effect.}
#'
#'   \item{get.tree}{ A single boolean value specifying whether to include
#'   available phylogenetic trees in the \code{TreeSummarizedExperiment}
#'   object. (By default: \code{get.tree = TRUE})}
#'
#'   \item{as.df}{ A single boolean value enabled when \code{output = "list"}.
#'   The argument specifies whether return functional data as a named list
#'   (one entry per element in the output list) of data.frames,
#'   with each data.frame containing results for all requested accessions.
#'   If \code{FALSE}, The function returns a list of lists, each element
#'   consisting of results for a single accession.
#'   (By default: \code{as.df = TRUE})}
#'
#'   \item{bulk.dl}{ A single boolean value specifying should
#'   MGnifyR attempt to speed things up by downloading
#'   relevant studies TSV results and only extracting the required columns,
#'   rather than using the JSONAPI interface. When getting results where
#'   multiple accessions share the same study, this option may result in
#'   significantly faster processing. However, there appear to be (quite a few)
#'   cases in the database where the TSV result columns do NOT match the
#'   expected accession names. This will hopefully be fixed in the future,
#'   but for now \code{bulk.dl} defaults to TRUE. When it does work, it can
#'   be orders of magnitude more efficient.
#'   (By default: \code{buld_dl = TRUE})}
#'
#' }
#'
#' @return
#' If only metagenomic data is retrieved, the result is returned in
#' \code{TreeSummarizedExperiment} object by default. The result can also be
#' returned as a \code{phyloseq} object or as a list of \code{data.frames}.
#' Note that \code{phyloseq} object can include only one phylogenetic tree
#' meaning that some taxa might be lost when data is subsetted based on tree.
#'
#' When functional data is retrieved in addition to metagenomic data, the result
#' is returned as a \code{MultiAssayExperiment} object. Other options are a list
#' containing \code{phyloseq} object and \code{data.frames} or just
#' \code{data.frames}.
#'
#' Functional data can be returned as a \code{MultiAssayExperiment} object or
#' as a list of \code{data.frames}.
#'
#' @examples
#' \dontrun{
#' # Get OTU tables as TreeSE
#' tse <- getResults(mg, accession_list, get.func=FALSE, get.taxa=TRUE)
#'
#' # Get functional data along with OTU tables as MAE
#' mae <- getResults(mg, accession_list, get.func=TRUE, get.taxa=TRUE)
#'
#' # Get same data as list
#' analyses_res <- getResults(
#'     mg, accession_list, get.func=TRUE, get.taxa=TRUE, output = "list",
#'     as.df=T, use.cache = T, bulk.dl = F)
#' }
#'
#' @name getResults
NULL

#' @rdname getResults
#' @include MgnifyClient.R utils.R
#' @importFrom plyr llply
#' @importFrom dplyr bind_rows
#' @importFrom reshape2 dcast
#' @importFrom stats as.formula
#' @export
setGeneric("getResults", signature = c("x"), function(
        x, accession, get.taxa = TRUE, get.func = TRUE,
        output = "TreeSE", use.cache = TRUE,
        ...
        )
    standardGeneric("getResults"))

#' @rdname getResults
#' @export
setMethod("getResults", signature = c(x = "MgnifyClient"), function(
        x, accession, get.taxa = TRUE, get.func = TRUE,
        output = "TreeSE", use.cache = TRUE,
        ...
        ){
    ############################### INPUT CHECK ################################
    if( !(.is_non_empty_character(accession)) ){
        stop("'accession' must be a single character value or list of ",
             "character values specifying the MGnify accession identifier.",
             call. = FALSE)
    }
    # If only one value, create a vector from it
    if( length(accession) == 1 ){
        accession <- c(accession)
    }
    if( !.is_a_bool(get.taxa) ){
        stop("'get.taxa' must be TRUE or FALSE.",
             call. = FALSE)
    }
    if( !(.is_a_bool(get.func) ||
          (is.character(get.func) &&
           all(get.func %in% names(.analyses_results_type_parsers)))) ){
        stop("'get.func' must be TRUE or FALSE or a single character value ",
             "or a list of character values specifying functional analysis ",
             "types.", call. = FALSE)
    }
    # Get all values, if TRUE
    if(!is.character(get.func) && get.func){
        get.func <- names(.analyses_results_type_parsers)
    }
    if( !(length(output) == 1 && output %in% c("list", "phyloseq", "TreeSE")) ){
        stop("'output' must be a 'TreeSE', 'list' or 'phyloseq'.",
             call. = FALSE)
    }
    if( !.is_a_bool(use.cache) ){
        stop("'use.cache' must be a single boolean value specifying whether to ",
             "use on-disk caching.", call. = FALSE)
    }
    ############################# INPUT CHECK END ##############################
    # Get functional data if user specified
    if( is.character(get.func) ){
        # If single value, create a list from it
        if( length(get.func) == 1){
            get.func <- c(get.func)
        }
        func_res <- .mgnify_get_analyses_results(
            client = x, accession = accession, retrievelist = get.func,
            output = output, use.cache = use.cache, ...
        )
    } else{
        func_res <- NULL
    }
    # Get microbial profiling data
    if( get.taxa ){
        taxa_res <- .mgnify_get_analyses_treese(
            client = x, accession = accession, use.cache = use.cache, ...
        )
    } else{
        taxa_res <- NULL
    }
    # Convert results to specified output type
    if( output != "list" ){
        result <- .convert_results_to_object(
            taxa_res, func_res, output, accession = accession)
    } else{
        # Create a final result list, if output is specified to be a list
        result <- append(taxa_res, func_res)
    }
    return(result)
})

################################ HELP FUNCTIONS ################################

# Convert results to TreeSE, MultiAssayExperiment or phyloseq.
#' @importFrom methods is
.convert_results_to_object <- function(taxa_res, func_res, output, accession){
    # If there are microbial profiling data, convert it to TreeSE or phyloseq
    if( !is.null(taxa_res) ){
        # Get sample metadata
        col_data <- taxa_res$sample_metadata
        # Bind sample metadata to one table
        col_data <- do.call(rbind, col_data)
        col_data <- DataFrame(col_data)
        # Merge TreeSEs into one
        tse_list <- taxa_res$tse_objects
        # Merge individual TreeSEs into one
        result <- mergeSEs(tse_list, assay.type = "counts", missing_values = 0)
        # Replace sample names with sample metadata sample names
        # (They are correct since the data was fetched in same order as TreeSEs)
        colnames(result) <- rownames(col_data)
        # Add sample metadata to the object
        colData(result) <- col_data
        # If user  wants phyloseq, convert TreeSE
        if( output == "phyloseq" ){
            result <- makePhyloseqFromTreeSE(result)
        }

    } else{
        result <- NULL
    }
    # If there are functional data
    if( !is.null(func_res) ){
        # If user wants TreeSE
        if( output == "TreeSE" ){
            # Create TreeSE from functional data
            func_res <- lapply(
                func_res, .create_TreeSE_from_func_data,
                accession = accession)
            # Remove data that is NULL
            func_res <- func_res[!unlist(lapply(func_res, is.null))]
            # Get colData from microbial profiling data TreeSE,
            # if it is included
            args <- list()
            if( !is.null(result) ){
                col_data <- colData(result)
                # Create a MAE
                result <- list(microbiota = result)
                exp_list <- append(result, func_res)
                args$colData <- col_data
            } else {
                exp_list <- func_res
                col_data <- NULL
            }
            exp_list <- ExperimentList(exp_list)
            args$experiments <- exp_list
            result <- do.call(MultiAssayExperiment, args)
        } else{
            # If user wants output as a phyloseq, give a list of one phyloseq
            # object and functional data
            result <- list(microbiota = result)
            result <- append(result, func_res)
        }
    }
    return(result)
}

# Create a TreeSE from single functional data data.frame
#' @importFrom S4Vectors SimpleList
.create_TreeSE_from_func_data <- function(x, accession){
    # If data was provided
    if( !is.null(x) ){
        # Add rownames
        rownames(x) <- x$accession
        # Get assay
        assay <- x[, colnames(x) %in% accession, drop = FALSE]
        # Get row_data
        row_data <- x[, !colnames(x) %in% accession]
        # Create a TreeSE
        assay <- as.matrix(assay)
        assays <- SimpleList(counts = assay)
        row_data <- DataFrame(row_data)
        tse <- TreeSummarizedExperiment(assays = assays, rowData = row_data)
    } else{
        tse <- NULL
    }
    return(tse)
}

# Helper function for importing microbial profiling data.
.mgnify_get_analyses_treese <- function(
        client, accession, use.cache,
        taxa.su = "SSU", get.tree = FALSE, ...){
    ############################### INPUT CHECK ################################
    if( !(.is_non_empty_string(taxa.su)) ){
        stop("'taxa.su' must be a single character value specifying taxa ",
             "subunit.",
             call. = FALSE)
    }
    if( !.is_a_bool(get.tree) ){
        stop("'get.tree' must be TRUE or FALSE.",
             call. = FALSE)
    }
    ############################# INPUT CHECK END ##############################
    # Some biom files don't import - so need a try/catch
    tse_list <- llply(accession, function(x) {
        tryCatch(
            .mgnify_get_single_analysis_treese(
                client, x, use.cache = use.cache, taxa.su = taxa.su,
                get.tree = get.tree),
            error=function(x){
                message("a biom listed in \"accession\" is missing from the ",
                        "retrieved tse_list")
                NULL}
        )
    }, .progress = "text")
    # The sample_data has been corrupted by doing the merge (names get messed
    # up and duplicated), so just regrab it with another lapply/rbind
    col_data <- lapply(accession, function(x){
        .mgnify_get_single_analysis_metadata(client, x, use.cache = use.cache)})
    # If user wants result as list
    result <- list(tse_objects=tse_list, sample_metadata = col_data)
    return(result)
}

################################ HELP FUNCTIONS ################################

# Get a single biom file and convert it to TreeSummarizedExperiment format
#
#' @importFrom mia loadFromBiom
#' @importFrom mia checkTaxonomy
#' @importFrom urltools parameters parameters<-
#' @importFrom httr GET
#' @importFrom httr write_disk
#' @importFrom ape read.tree
#' @importFrom TreeSummarizedExperiment rowTree
.mgnify_get_single_analysis_treese <- function(
        client = NULL, accession, use.cache = TRUE, downloadDIR = NULL,
        taxa.su = "SSU", get.tree = FALSE){
    # Get the metadata
    metadata_df <- .mgnify_get_single_analysis_metadata(
        client, accession, use.cache=use.cache)
    analysis_data <- .mgnify_retrieve_json(
        client, paste("analyses",accession, sep="/"), use.cache = use.cache)
    download_url <- analysis_data[[1]]$relationships$downloads$links$related
    analysis_downloads <- .mgnify_retrieve_json(
        client, complete_url = download_url, use.cache = use.cache)

    # Depending on the pipeline version, there may be more than one OTU table
    # available (LSU/SSU), so try and get the one specified in taxa.su -
    # otherwise spit out a warning and grab the generic (older pipelines)
    available_biom_files <- analysis_downloads[grepl('JSON Biom', sapply(
        analysis_downloads, function(x){x$attributes$`file-format`$name}))]
    biom_position <- grepl(taxa.su, sapply(
        available_biom_files, function(x){x$attributes$`group-type`}))
    if(sum(biom_position) == 0){
        if(client@warnings){
            warning("Unable to locate requested taxonomy type ", taxa.su, ". ",
                    "This is likely due to the current analysis having been ",
                    "performed on an older version of the MGnify pipeline. ",
                    "The available BIOM file will be used instead.")
        }
        biom_url <- available_biom_files[[1]]$links$self
    }else{
        biom_url <- available_biom_files[biom_position][[1]]$links$self
    }

    # Can specify a seperate dir for saving biom files, otherwise they end up
    # in the client@cachdir folder, under "bioms"
    if (is.null(downloadDIR)){
        downloadDIR <- paste(client@cacheDir,"biom_files",sep="/")
        dir.create(downloadDIR, recursive = T, showWarnings = client@warnings)
    }
    # Clear out any ?params after the main location - don't need them for this
    parameters(biom_url) <- NULL

    fname <- utils::tail(strsplit(biom_url, '/')[[1]], n=1)
    biom_path <- paste(downloadDIR, fname, sep="/")

    ## Quick check to see if we should clear the disk cache ~for this specific
    # call~ - used for debugging and when MGnify breaks
    if(use.cache && client@clearCache){
        message(paste("clear_cache is TRUE: deleting ", biom_path, sep=""))
        tryCatch(unlink(biom_path), error=warning)
    }

    if (! file.exists(biom_path)){#} | !use_downloads ){
        GET(biom_url, write_disk(biom_path, overwrite = TRUE))
    }
    # Load in the TreeSummarizedExperiment object
    # Some files do not have sample names --> suppress warning message
    # "there is no colnames, you can add them..."
    suppressWarnings(
        tse <- loadFromBiom(biom_path, removeTaxaPrefixes = TRUE,
                            rankFromPrefix = TRUE
                            )
    )
    # Add sample names if data does not include them
    if( is.null(colnames(tse)) ){
        colnames(tse) <- paste0("sample_", seq_len(ncol(tse)))
    }

    # # Need to check if the taxonomy was parsed correctly - depending on the
    # # pipeline it may need a bit of help:
    # checkTaxonomy(tse)
    # # if (ncol(tax_table(psobj)) == 1){
    # #     psobj <- import_biom(biom_path, parseFunction = parse_taxonomy_qiime)
    # # }
    # # if(! "Kingdom" %in% colnames(tax_table(psobj))){
    # #     psobj <- import_biom(
    # #         biom_path, parseFunction = parse_taxonomy_greengenes)
    # # }
    if(get.tree){
        # Is there a tree?
        tvec <- grepl('Phylogenetic tree', sapply(
            analysis_downloads, function(x) x$attributes$`description`$label))
        if(any(tvec)){
            tree_url <- analysis_downloads[tvec][[1]]$links$self
            #Clear out any ?params after the main location - don't need them for this
            parameters(tree_url) <- NULL

            fname <- utils::tail(strsplit(tree_url, '/')[[1]], n=1)
            tree_path <- paste(downloadDIR, fname, sep="/")

            ## Quick check to see if we should clear the disk cache
            # ~for this specific call~ - used for debugging and when MGnify breaks
            if(use.cache && client@clearCache){
                message(paste("clear_cache is TRUE: deleting ",tree_path, sep=""))
                tryCatch({unlink(tree_path)}, error=warning)
            }

            if (! file.exists(tree_path)){#} | !use_downloads ){
                GET(tree_url, write_disk(tree_path, overwrite = T ))
            }
        }
        row_tree <- read.tree(tree_path)
        rowTree(tse) <- row_tree
    }
    return(tse)
}

# Helper function for importing functional data
.mgnify_get_analyses_results <- function(
        client, accession, retrievelist, output, use.cache,
        as.df = TRUE, bulk.dl = TRUE, ...){
    ############################### INPUT CHECK ################################
    # If output is TreeSE or object, get results as data.frames
    as.df <- ifelse(output != "list", TRUE, as.df)
    if( !.is_a_bool(as.df) ){
        stop("'as.df' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(bulk.dl) ){
        stop("'bulk.dl' must be TRUE or FALSE.", call. = FALSE)
    }
    ############################# INPUT CHECK END ##############################
    # Get functional results
    all_results <- llply(accession, .progress = "text", function(x){
        .mgnify_get_single_analysis_results(
            client, x, use.cache = use.cache, retrievelist = retrievelist,
            bulk.files = bulk.dl)
    })
    # Add names based on accessions codes
    names(all_results) <- accession

    # Compact the result type dataframes into a single instance. Per
    # accession counts in each column.
    if(as.df){
        # Create a list data.frames from results
        all_results <- llply(retrievelist, function(y){
            tryCatch({
                r <- lapply(all_results, function(x){
                    df <- as.data.frame(x[[y]])
                    df
                })
                longform <- bind_rows(r, .id = "analysis")
                cn <- colnames(longform)
                extras <- cn[!(cn %in% c("count","index_id", "analysis"))]
                final_df <- dcast(
                    longform, as.formula(paste(paste(extras, collapse = " + "),
                                               " ~ analysis")),
                    value.var = "count", fun.aggregate = sum)
                final_df}, error=function(x) NULL)
        })
        # Give names based on functional analysis types
        names(all_results) <- retrievelist
    }
    return(all_results)
}

# types can be found using names(MGnifyR:::.analyses_results_type_parsers).
# Note that not depending on the particular analysis type, pipeline
# Retrieves combined study/sample/analysis metadata - not exported
#' @importFrom urltools parameters parameters<-
#' @importFrom httr GET
#' @importFrom httr write_disk
#' @importFrom dplyr bind_rows
#' @importFrom utils read.csv2
.mgnify_get_single_analysis_results <- function(
        client = NULL, accession, retrievelist=c(), use.cache = TRUE,
        max.hits = -1, bulk.files = FALSE){
    metadata_df <- .mgnify_get_single_analysis_metadata(
        client, accession, use.cache=use.cache, max.hits = max.hits)

    # Should we try and grab the study's full TSV download rather than parse
    # through the JSON API? Doing so has the potential to use a LOT more disk
    # space, along with potentially increased data download. It should
    # be faster though, except in pathological cases (e.g. only 1 sample per
    # 1000 sample study required). As with everything else, we make use of
    # local caching to speed things along.
    if(bulk.files){
        downloadDIR <- paste(client@cacheDir, "tsv", sep="/")
        if(!dir.exists(downloadDIR)){
            dir.create(
                downloadDIR, recursive = TRUE, showWarnings = client@warnings)
        }

        available_downloads_json <- .mgnify_retrieve_json(
            client, path=paste("studies", metadata_df$study_accession,
                               "downloads", sep="/"), use.cache = use.cache)
        #cat(str(dput(available_downloads_json)))
        parsed_results <- lapply(available_downloads_json, function(r) {
            # Figure out the label mapping
            cur_lab <- r$attributes$description$label

            cur_pipeversion <- r$relationships$pipeline$data$id

            # There MUST be a better way to do the line below...
            cur_type <- names(
                .analyses_results_bulk_file_names[is.finite(match(
                    .analyses_results_bulk_file_names, cur_lab))])

            if (!identical(cur_type, character(0))){

                # Check the pipeline versions match
                if (cur_pipeversion == metadata_df$`analysis_pipeline-version`){
                    if( cur_type %in% retrievelist) {

                        #Get the url
                        data_url <- r$links$self

                        #Clear off extraneous gubbins
                        parameters(data_url) <- NULL

                        #build the cache filename
                        fname <- utils::tail(strsplit(data_url, '/')[[1]], n=1)

                        #At this point we might have alread got the data we want
                        # loaded. Check the memory cache object

                        if(client@useMemCache & (
                            cur_type %in% names(mgnify_memory_cache)) & (
                                mgnify_memory_cache[cur_type][
                                    "fname"] == fname)){
                            tmp_df <- mgnify_memory_cache[cur_type][["data"]]
                        }else{
                            # Nope - gonna have to load it up from disk or grab
                            # it from t'interweb
                            data_path <- paste(downloadDIR, fname, sep="/")

                            if(use.cache & client@clearCache){
                                message(paste("clear_cache is TRUE: deleting ",
                                              data_path, sep=""))
                                tryCatch({unlink(data_path)}, error=warning)
                            }

                            if (! file.exists(data_path)){#} | !use_downloads ){
                                GET(data_url, write_disk(
                                    data_path, overwrite = TRUE ))
                            }

                            # Load the file (might be big so save it in the
                            # 1 deep cache)
                            tmp_df <- read.csv2(
                                data_path, sep="\t", header = TRUE,
                                stringsAsFactors = FALSE)
                        }
                        # Save it in memory using "super assignment" - which
                        # I'm not really sure about but it seems to work...
                        # thing'd be much easier if R passed objects
                        # by reference.
                        if(client@useMemCache){
                            mgnify_memory_cache[[cur_type]] <- list(
                                data=tmp_df, fname=fname)
                        }

                        # Because there seem to be "mismatches" between the
                        # JSON and downloadable files, and maybe some issues
                        # with missing downloads, we have to check if we
                        # actually got a valid file:
                        if(ncol(tmp_df) < 3){
                            warning(paste(
                                "Invalid download for", accession, sep=" "))
                            return(NULL)
                        }

                        # Need to figure out how many columns to keep -
                        # the first one is always an ID, but need to keep some
                        # others as well...
                        i <- 1
                        #tmp_df <- read.csv('~/.MGnify_cache/tsv/ERP108138_
                        # IPR_abundances_v4.1.tsv', sep="\t", header = T,
                        # stringsAsFactors = F)
                        while(
                            any(is.na(suppressWarnings(
                                as.numeric(tmp_df[,i]))))
                            ){
                            i <- i+1
                        }
                        i <- i-1

                        # Also need the column name for this particular
                        # analysis...
                        # As far as I can see they could be either assembly IDs
                        # or run ids. FFS. Assuming that both assembly and run
                        # won't be present...:
                        if("assembly_accession" %in% colnames(metadata_df)){
                            accession <- metadata_df$assembly_accession[[1]]
                        }else if("run_accession" %in% colnames(metadata_df)){
                            accession <- metadata_df$run_accession[[1]]
                        }
                        #cat(accession)
                        #    cat(str(head(tmp_df[1:5,1:5])))
                        # Break up the dataframe, only keeping the bits we need.

                        # So at this point we learn that some of the "download"
                        # files don't match the assembly/run IDs given in the
                        # JSON. For now, do a try/catch, chuck a warning, and
                        # then optionally go off and try again - this time from
                        # the JSON. No doubt this'll be fixed at some point in
                        # the future...

                        column_position <- match(accession, colnames(tmp_df))
                        if (is.na(column_position)){
                            warning(paste(
                                "Failed to find column", accession, sep = " "))
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
        # R is sometimes a bit ~awkward~
        names(parsed_results) <- names(
            .analyses_results_bulk_file_names)[match(unlist(lapply(
                available_downloads_json,
                function(x){x$attributes$description$label})),
                .analyses_results_bulk_file_names)]
    }else{
        # Now (re)load the analysis data:
        analysis_data <- .mgnify_retrieve_json(
            client, paste("analyses", accession, sep="/"),
            use.cache = use.cache, max.hits = max.hits)
        # For now try and grab them all - just return the list -
        # don't do any processing...
        all_results <- lapply(names(
            .analyses_results_type_parsers), function(r) {
            if(r %in% retrievelist){
                tmp <- .mgnify_retrieve_json(
                    client,
                    complete_url = analysis_data[[1]]$relationships[[
                        r]]$links$related,
                    use.cache = use.cache, max.hits = max.hits)
                tmp
            }
        })
        names(all_results) <- names(.analyses_results_type_parsers)
        parsed_results <- sapply(names(all_results), function(x){
            all_json <- all_results[[x]]
            if(! is.null(all_json)){
                res_df <- do.call(
                    bind_rows, lapply(
                        all_json,.analyses_results_type_parsers[[x]] ))
                rownames(res_df) <- res_df$index_id
                res_df
            }else{
                NULL
            }
        })
    }
    # Return the results...
    parsed_results
}

# Result table caching
mgnify_memory_cache <- list()

# Which parser do you use for which type of output?
# a list of parsers for each output type.
.analyses_results_type_parsers <- list(
    `taxonomy` = .mgnify_parse_tax,
    `taxonomy-itsonedb` = .mgnify_parse_tax,
    `go-slim`=.mgnify_parse_func,
    `taxonomy-itsunite` = .mgnify_parse_tax,
    `taxonomy-ssu` = .mgnify_parse_tax,
    `taxonomy-lsu` = .mgnify_parse_tax,
    `antismash-gene-clusters` = .mgnify_parse_func,
    `go-terms` = .mgnify_parse_func,
    `interpro-identifiers` = .mgnify_parse_func)

# This maps the json attribute name for retrievelist to the "description ->
# label" attribute in the study downloads section
.analyses_results_bulk_file_names <- list(
    `taxonomy` = "Taxonomic assignments SSU",
    `taxonomy-itsonedb` = "Taxonomic assignments ITS",
    `go-slim` = "GO slim annotation",
    `taxonomy-itsunite` = "Taxonomic assignments Unite",
    `taxonomy-ssu` = "Taxonomic assignments SSU",
    `taxonomy-lsu` = "Taxonomic assignments LSU",
    `antismash-gene-clusters` = .mgnify_parse_func,
    `go-terms` = "Complete GO annotation",
    `interpro-identifiers` = "InterPro matches",
    `phylo-tax-ssu` = "Phylum level taxonomies SSU",
    `phylo-tax-lsu` = "Phylum level taxonomies LSU")
