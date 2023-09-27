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
#' functional results will be available. (By default: \code{get.func = TRUE})
#'
#' @param use.cache A single boolean value specifying whether to use the
#' MGnify local caching system to speed up searching. It is highly
#' recommended that this be enabled. Note that files are downloaded to local system
#' when they are fetched from the database. The files are not removed meaning
#' that the local storage can include additional files after the run even though
#' \code{use.cache = FALSE} was specified. (By default: \code{use.cache = TRUE})
#'
#' @param verbose A single boolean value to specify whether to show
#' the progress bar. (By default: \code{verbose = TRUE})
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
#' tse <- getResult(mg, accession_list, get.func=FALSE, get.taxa=TRUE)
#'
#' # Get functional data along with OTU tables as MAE
#' mae <- getResult(mg, accession_list, get.func=TRUE, get.taxa=TRUE)
#'
#' # Get same data as list
#' analyses_res <- getResult(
#'     mg, accession_list, get.func=TRUE, get.taxa=TRUE, output = "list",
#'     as.df=T, use.cache = T, bulk.dl = F)
#' }
#'
#' @name getResult
NULL

#' @rdname getResult
#' @include MgnifyClient.R utils.R
#' @importFrom plyr llply
#' @importFrom dplyr bind_rows
#' @importFrom reshape2 dcast
#' @importFrom stats as.formula
#' @export
setGeneric("getResult", signature = c("x"), function(
        x, accession, get.taxa = TRUE, get.func = TRUE,
        output = "TreeSE", use.cache = TRUE, verbose = TRUE,
        ...
        )
    standardGeneric("getResult"))

#' @rdname getResult
#' @export
setMethod("getResult", signature = c(x = "MgnifyClient"), function(
        x, accession, get.taxa = TRUE, get.func = TRUE,
        output = "TreeSE", use.cache = TRUE, verbose = TRUE,
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
    if( !.is_a_bool(verbose) ){
        stop("'verbose' must be a single boolean value specifying whether to ",
             "show progress.", call. = FALSE)
    }
    verbose <- ifelse(verbose, "text", "none")
    ############################# INPUT CHECK END ##############################
    # Get functional data if user specified
    if( is.character(get.func) ){
        # If single value, create a vector from it
        if( length(get.func) == 1){
            get.func <- c(get.func)
        }
        func_res <- .mgnify_get_analyses_results(
            client = x, accession = accession, retrievelist = get.func,
            output = output, use.cache = use.cache, verbose = verbose, ...
        )
    } else{
        func_res <- NULL
    }
    # Get microbial profiling data
    if( get.taxa ){
        taxa_res <- .mgnify_get_analyses_treese(
            client = x, accession = accession, use.cache = use.cache,
            verbose = verbose, ...
        )
    } else{
        taxa_res <- NULL
    }
    # Convert results to specified output type
    if( output != "list" ){
        result <- .convert_results_to_object(
            taxa_res, func_res, output)
    } else{
        # Create a final result list, if output is specified to be a list
        result <- append(taxa_res, func_res)
    }
    return(result)
})

################################ HELP FUNCTIONS ################################

# Convert results to TreeSE, MultiAssayExperiment or phyloseq.
#' @importFrom methods is
#' @importFrom dplyr bind_rows
.convert_results_to_object <- function(taxa_res, func_res, output){
    result <- NULL
    # If there are microbial profiling data, convert it to TreeSE or phyloseq
    if( !is.null(taxa_res) ){
        # Get TreeSE objects
        tse_list <- taxa_res$tse_objects
        # Get sample metadata
        col_data <- taxa_res$sample_metadata
        # If some results were not found, remove them
        ind <- !unlist(lapply(tse_list, is.null))
        # If there are samples left after subsetting
        if( any(ind) ){
            tse_list <- tse_list[ind]
            col_data <- col_data[ind]
            # Bind sample metadata to one table
            col_data <- do.call(bind_rows, col_data)
            col_data <- DataFrame(col_data)
            # Merge individual TreeSEs into one
            result <- mergeSEs(tse_list, assay.type = "counts", missing_values = 0)
            # Order the sample metadata
            col_data <- col_data[ colnames(result), , drop = FALSE]
            # Add sample metadata to the object
            colData(result) <- col_data
            # If user wants phyloseq, convert TreeSE
            if( output == "phyloseq" ){
                result <- makePhyloseqFromTreeSE(result)
            }
        } else{
            warning("No taxonomy data was found for the dataset.", call. = FALSE)
        }
    }
    # If there are functional data
    if( !is.null(func_res) ){
        # Remove data that is NULL
        ind <- !unlist(lapply(func_res, is.null))
        # If there are experiments left after subsetting
        if( any(ind) ){
            # Take only those experiments that are not NULL
            func_res <- func_res[ind]
            # If user wants TreeSE
            if( output == "TreeSE" ){
                # Create TreeSE from functional data
                func_res <- lapply(
                    func_res,
                    .create_TreeSE_from_func_data, tse = result)
                # Get colData from microbial profiling data TreeSE,
                # if it is included
                args <- list()
                # If taxa data is included
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

                # If there are more than 1 experiments, create MAE
                if( length(exp_list) > 1 ){
                    exp_list <- ExperimentList(exp_list)
                    args$experiments <- exp_list
                    result <- do.call(MultiAssayExperiment, args)
                } else{
                    # If there are only 1 experiment, give it as it is
                    result <- exp_list[[1]]
                }
            } else{
                # If user wants output as a phyloseq, give a list of one phyloseq
                # object and functional data
                result <- list(microbiota = result)
                result <- append(result, func_res)
                # If there are only one experiment, take it out from the list
                if( length(result) == 1 ){
                    result <- result[[1]]
                }
            }
        } else{
            warning("No functional data found for the dataset.", call. = FALSE)
        }
    }
    return(result)
}

# Create a TreeSE from single functional data data.frame
#' @importFrom S4Vectors SimpleList
#' @importFrom dplyr %>% mutate_all na_if
.create_TreeSE_from_func_data <- function(x, tse_microbiota){
    # If data was provided
    if( !is.null(x) ){
        # Get assay
        assay <- dcast(x, index_id ~ analysis, value.var = "count")
        rownames(assay) <- assay[["index_id"]]
        assay[["index_id"]] <- NULL
        # Get row_data
        row_data <- x[ , !colnames(x) %in% c("analysis", "count"), drop = FALSE]
        row_data <- row_data[!duplicated(row_data), , drop = FALSE]
        rownames(row_data) <- row_data[["index_id"]]
        row_data <- row_data[rownames(assay), , drop=FALSE]
        # Get taxonomy from the information (e.g. when the data is taxonomy)
        tax_tab <- mia:::.parse_taxonomy(row_data, column_name = "index_id")
        # If taxonomy information was found
        if( ncol(tax_tab) > 0 ){
            # Remove prefixes
            tax_tab <- mia:::.remove_prefixes_from_taxa(tax_tab)
            # Replace empty cells with NA
            tax_tab <- tax_tab %>% as.data.frame() %>% mutate_all(na_if, "")
            # Add taxonomy info to original rowData
            row_data <- cbind(row_data, tax_tab)
        }
        # If microbiota data exists, order the functional data based on that
        # Do not drop samples that are not found from microbiota data
        if( !is.null(tse_microbiota) ){
            assay <-  assay[ , order(match(
                colnames(assay), colnames(tse_microbiota))), drop = FALSE]
        }
        # Create a TreeSE
        assay <- as.matrix(assay)
        assays <- SimpleList(counts = assay)
        row_data <- DataFrame(row_data)
        # Get arguments for TreeSE
        args <- list(assays = assays, rowData = row_data)
        if( !is.null(tse_microbiota) ){
            args$colData <- colData(tse_microbiota)
        }
        # Create TreeSE
        tse <- do.call(TreeSummarizedExperiment, args)
    } else{
        tse <- NULL
    }
    return(tse)
}

# Helper function for importing microbial profiling data.
.mgnify_get_analyses_treese <- function(
        client, accession, use.cache, verbose,
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
    # Give message about progress
    if( verbose =="text" ){
        message("Fetching taxonomy data...")
    }
    # Get TreeSE objects
    tse_list <- llply(accession, function(x) {
            .mgnify_get_single_analysis_treese(
                client, x, use.cache = use.cache, taxa.su = taxa.su,
                get.tree = get.tree, ...)
    }, .progress = verbose)
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
        taxa.su = "SSU", get.tree = FALSE, ...){
    # Get the metadata related to analysis
    metadata_df <- .mgnify_get_single_analysis_metadata(
        client, accession, use.cache=use.cache, ...)
    analysis_data <- .mgnify_retrieve_json(
        client, paste("analyses",accession, sep="/"), use.cache = use.cache, ...)
    # Get the metadata related to available files etc
    download_url <- analysis_data[[1]]$relationships$downloads$links$related
    analysis_downloads <- .mgnify_retrieve_json(
        client, complete_url = download_url, use.cache = use.cache, ...)

    # Depending on the pipeline version, there may be more than one OTU table
    # available (LSU/SSU), so try and get the one specified in taxa.su -
    # otherwise spit out a warning and grab the generic (older pipelines)
    available_biom_files <- analysis_downloads[grepl('JSON Biom', sapply(
        analysis_downloads, function(x){x$attributes$`file-format`$name}))]
    # Check if any biom files was found
    if( is.null(available_biom_files) || length(available_biom_files) == 0 ){
        warning("No BIOM data found for accession '", accession, "'.", call. = FALSE)
        return(NULL)
    }
    biom_position <- grepl(taxa.su, sapply(
        available_biom_files, function(x){x$attributes$`group-type`}))
    if( sum(biom_position) == 0 ){
        if( client@warnings ){
            warning("Unable to locate requested taxonomy type ", taxa.su, ". ",
                    "This is likely due to the current analysis having been ",
                    "performed on an older version of the MGnify pipeline. ",
                    "The available BIOM file will be used instead.",
                    call. = FALSE)
        }
        biom_url <- available_biom_files[[1]]$links$self
    } else {
        biom_url <- available_biom_files[biom_position][[1]]$links$self
    }

    # Can specify a separate dir for saving biom files, otherwise they end up
    # in the client@cachdir folder, under "bioms"
    if(is.null(downloadDIR)){
        downloadDIR <- paste(client@cacheDir,"biom_files",sep="/")
        dir.create(downloadDIR, recursive = T, showWarnings = client@warnings)
    }
    # Clear out any ?params after the main location - don't need them for this
    parameters(biom_url) <- NULL
    # Get the file name and path to it in local machine
    fname <- utils::tail(strsplit(biom_url, '/')[[1]], n=1)
    biom_path <- paste(downloadDIR, fname, sep="/")

    ## Quick check to see if we should clear the disk cache  for this specific
    # call  - used for debugging and when MGnify breaks
    if( use.cache && client@clearCache && file.exists(biom_path) ){
        message(paste("clear_cache is TRUE: deleting ", biom_path, sep=""))
        unlink(biom_path)
    }
    # Download the file from the database to specific file path
    fetched_from_url <- FALSE
    if( !file.exists(biom_path) ){
        res <- GET(biom_url, write_disk(biom_path, overwrite = TRUE))
        fetched_from_url <- TRUE
        # If the file was not successfully downloaded
        if( res$status_code != 200 ){
            warning(
                biom_url, ": ", content(res, ...)$errors[[1]]$detail,
                " A biom listed in 'accession' is missing from the ",
                "output.", call. = FALSE)
            # Remove the downloaded file, it includes only info on errors
            unlink(biom_path)
            return(NULL)
        }
    }

    # Load in the TreeSummarizedExperiment object
    tse <- loadFromBiom(
        biom_path, removeTaxaPrefixes = TRUE, rankFromPrefix = TRUE,
        remove.artifacts = TRUE)
    # If the file was not in store already but fetched from database, and cache
    # storing is disabled
    if( fetched_from_url && !use.cache ){
        unlink(biom_path)
    }
    # TreeSE has sample ID as its colnames. Rename so that it is the accession ID.
    colData(tse)[["biom_sample_id"]] <- colnames(tse)
    colnames(tse) <- accession

    # If user wants also phylogenetic tree
    if(get.tree){
        # Is there a tree?
        tvec <- grepl('Phylogenetic tree', sapply(
            analysis_downloads, function(x) x$attributes$`description`$label))
        if( any(tvec) ){
            # Get the url address of tree
            tree_url <- analysis_downloads[tvec][[1]]$links$self
            # Clear out any ?params after the main location - don't need them for this
            parameters(tree_url) <- NULL
            # Get the file name for the tree
            fname <- utils::tail(strsplit(tree_url, '/')[[1]], n=1)
            # Get the path for the tree
            tree_path <- paste(downloadDIR, fname, sep="/")

            ## Quick check to see if we should clear the disk cache
            #  for this specific call  - used for debugging and when MGnify breaks
            if(use.cache && client@clearCache && file.exists(tree_path) ){
                message(paste("clear_cache is TRUE: deleting ",
                              tree_path, sep=""))
                unlink(tree_path)
            }
            # Download the file from the database to specific file path
            fetched_from_url <- FALSE
            if( !file.exists(tree_path) ){
                res <- GET(tree_url, write_disk(tree_path, overwrite = TRUE))
                fetched_from_url <- TRUE
                # If the file was not successfully downloaded
                if( res$status_code != 200 ){
                    warning(
                        tree_url, ": ", content(res, ...)$errors[[1]]$detail,
                        " A phylogenetic tree listed in 'accession' is ",
                        "missing from the output.", call. = FALSE)
                }
            }
            # Add the tree to TreeSE object
            if( file.exists(tree_path) ){
                row_tree <- read.tree(tree_path)
                rowTree(tse) <- row_tree
                # If the file was not in store already but fetched from database,
                # and cache storing is disabled
                if( fetched_from_url && !use.cache ){
                    unlink(biom_path)
                }
            }
        }
    }
    return(tse)
}

# Helper function for importing functional data
.mgnify_get_analyses_results <- function(
        client, accession, retrievelist, output, use.cache, verbose,
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
    # Give message about progress
    if( verbose == "text" ){
        message("Fetching functional data...")
    }
    # Get functional results
    all_results <- llply(accession, function(x){
        .mgnify_get_single_analysis_results(
            client, x, use.cache = use.cache, retrievelist = retrievelist,
            bulk.files = bulk.dl, ...)
    }, .progress = verbose)
    # Add names based on accessions codes
    names(all_results) <- accession

    # Compact the result type dataframes into a single instance. Per
    # accession counts in each column.
    if(as.df){
        # Check which data types are available
        data_types <- names(all_results[[1]])
        # Combine results data type -wise
        all_results <- lapply(data_types, function(data_type){
            # Get only certain type of data of all samples
            temp <- lapply(all_results, function(acc_res) acc_res[[data_type]])
            # Combine data to df
            temp <- bind_rows(temp, .id = "analysis")
            # If there was no data, the df is empty
            if( nrow(temp) > 0 && ncol(temp) > 0 ){
                # Counts are numeric values, convert...
                temp[["count"]] <- as.numeric(temp[["count"]])
                # There are three columns that we want; "analysis"
                # (sample ID), "count" (counts), and "index_id" (feature)
                # columns.

                # If the "index_id" column equals with the "accession" column,
                # drop the latter. They should match in taxonomy data, but
                # functional might have additional columns.
                if( all(temp[["index_id"]] == temp[["accession"]]) ){
                    temp <- temp[ , colnames(temp) != "accession", drop = FALSE]
                }
            } else {
                temp <- NULL
            }
            return(temp)
        })
        # Give names based on data types
        names(all_results) <- data_types
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
        max.hits = NULL, bulk.files = FALSE, ...){
    # Get the metadata describing the samples
    metadata_df <- .mgnify_get_single_analysis_metadata(
        client, accession, use.cache=use.cache, max.hits = max.hits)

    # Should we try and grab the study's full TSV download rather than parse
    # through the JSON API? Doing so has the potential to use a LOT more disk
    # space, along with potentially increased data download. It should
    # be faster though, except in pathological cases (e.g. only 1 sample per
    # 1000 sample study required). As with everything else, we make use of
    # local caching to speed things along.
    if(bulk.files){
        # Create a directory for donwloading the file
        downloadDIR <- paste(client@cacheDir, "tsv", sep="/")
        if(!dir.exists(downloadDIR)){
            dir.create(
                downloadDIR, recursive = TRUE, showWarnings = client@warnings)
        }
        # Get what files are available in database
        available_downloads_json <- .mgnify_retrieve_json(
            client,
            path=paste("studies", metadata_df$study_accession, "downloads", sep="/"),
            use.cache = use.cache, ...)
        # Load the files
        parsed_results <- lapply(available_downloads_json, function(r) {
            # Figure out the label mapping
            cur_lab <- r$attributes$description$label
            # Pipe version
            cur_pipeversion <- r$relationships$pipeline$data$id

            # There MUST be a better way to do the line below...
            # Get the type of the data
            cur_type <- names(.analyses_results_bulk_file_names)[
                cur_lab == .analyses_results_bulk_file_names]

            # Check the pipeline versions match and load the file if the datatype
            # is specified to be loaded.
            if(length(cur_type) > 0 &&
               cur_pipeversion == metadata_df$`analysis_pipeline-version` &&
               any(cur_type %in% retrievelist)){
                # If there are "taxonomic assignments ssu", it can match with
                # 2 --> take the first one
                cur_type <- cur_type[[1]]
                # Get the data
                temp <- .get_bulk_files(
                    cur_type, client=client, r=r,
                    metadata_df=metadata_df, data_path=data_path,
                    downloadDIR=downloadDIR, use.cache=use.cache, ...)
            } else{
                temp <- NULL
            }
            return(list(type=cur_type, data=temp))
        })
        # Add data types to data as names (taxonomy might have 2 data types if NULL
        # because these both 2 data types are tried to fetch)
        cur_type <- unlist(lapply(parsed_results, function(x)
            ifelse( length(x$type) > 1, x$type[[1]], x$type)))
        parsed_results <- lapply(parsed_results, function(x) x$data)
        names(parsed_results) <- cur_type
    }else{
        # If user do not want to fetch bulk files
        # Now (re)load the analysis data:
        analysis_data <- .mgnify_retrieve_json(
            client, paste("analyses", accession, sep="/"),
            use.cache = use.cache, max.hits = max.hits, ...)
        # For now try and grab all data types that user has specified -
        # just return the list - don't do any processing...
        all_results <- lapply(
            names(.analyses_results_type_parsers), function(r) {
                if(r %in% retrievelist){
                    tmp <- .mgnify_retrieve_json(
                        client,
                        complete_url = analysis_data[[1]]$relationships[[r]]$links$related,
                        use.cache = use.cache,
                        max.hits = max.hits, ...)
                } else{
                    tmp <- NULL
                }
                return(tmp)
        })
        # Add data types as names
        names(all_results) <- names(.analyses_results_type_parsers)
        parsed_results <- sapply(names(all_results), function(x){
            # Get the specific type of data
            all_json <- all_results[[x]]
            # If specific type of data can be found
            if( !is.null(all_json) && length(all_json) > 0 ){
                # Get the specific type of data from all accessions
                all_json <- lapply(all_json, .analyses_results_type_parsers[[x]])
                # Combine
                res_df <- do.call(bind_rows, all_json)
                rownames(res_df) <- res_df$index_id
            }else{
                res_df <- NULL
            }
            return(res_df)
        })
    }
    return(parsed_results)
}

# Helper function for getting bulk file
.get_bulk_files <- function(
        cur_type, client, r, metadata_df, data_path, downloadDIR, use.cache, ...){
    # Get the url
    data_url <- r$links$self
    # Clear off extraneous gubbins
    parameters(data_url) <- NULL
    #build the cache filename
    fname <- utils::tail(strsplit(data_url, '/')[[1]], n=1)

    #At this point we might have alread got the data we want
    # loaded. Check the memory cache object

    if( (client@useMemCache &&
         cur_type %in% names(mgnify_memory_cache) &&
         mgnify_memory_cache[cur_type]["fname"] == fname) ){
        tmp_df <- mgnify_memory_cache[cur_type][["data"]]
    }else{
        # Nope - gonna have to load it up from disk or grab
        # it from t'interweb
        data_path <- paste(downloadDIR, fname, sep="/")
        # Clear cache if specified
        if(use.cache && client@clearCache && file.exists(data_path) ){
            message(paste("clear_cache is TRUE: deleting ",
                          data_path, sep=""))
            unlink(data_path)
        }
        # Download the file from the database to specific file path
        if( !file.exists(data_path) ){
            res <- GET(data_url, write_disk(data_path, overwrite = TRUE))
            # If the file was not successfully downloaded
            if( res$status_code != 200 ){
                warning(
                    data_url, ": ", content(res, ...)$errors[[1]]$detail,
                    " Error while loading the file from database. The data ",
                    "from the file is not included in the output.",
                    call. = FALSE)
                # Remove the downloaded file
                unlink(data_path)
                return(NULL)
            }
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
    if( ncol(tmp_df) < 3 ){
        warning(
            paste("Invalid download for", accession, sep=" "), call. = FALSE)
        return(NULL)
    }

    # Need to figure out how many columns to keep -
    # Get those columns that do not include numeric values but some info about
    # samples etc... First column includes always sample ID.
    info_cols <- apply(tmp_df, 2, function(x) is.na(suppressWarnings(as.numeric(x))) )
    info_cols <- colSums(info_cols) == nrow(info_cols)
    info_cols <- which(info_cols)
    info_cols <- unique(c(1, info_cols))

    # Also need the column name for this particular
    # analysis...
    # As far as I can see they could be either assembly IDs
    # or run ids. FFS. Assuming that both assembly and run
    # won't be present...:
    if( "assembly_accession" %in% colnames(metadata_df) ){
        accession <- metadata_df$assembly_accession[[1]]
    }else if( "run_accession" %in% colnames(metadata_df) ){
        accession <- metadata_df$run_accession[[1]]
    } else{
        warning(
            paste("Failed to data on", accession, sep = " "))
        return(NULL)
    }

    # Get the correct sample and subset the data so that it includes only
    # the sample and info columns
    column_position <- match(accession, colnames(tmp_df))
    if( is.na(column_position) || length(column_position) != 1 ){
        warning(
            paste("Failed to find column", accession, sep = " "), call. = FALSE)
        return(NULL)
    }
    keeper_columns <- c(info_cols, column_position)
    tmp_df <- tmp_df[, keeper_columns, drop = FALSE]
    # Adjust colnames rownames and add column
    colnames(tmp_df)[1] <- "accession"
    colnames(tmp_df)[ ncol(tmp_df) ] <- "count"
    tmp_df$index_id <- tmp_df$accession
    rownames(tmp_df) <- tmp_df$accession
    return(tmp_df)
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
