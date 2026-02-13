#' Get microbial and/or functional profiling data for a list of accessions
#'
#' @details
#' Given a set of analysis accessions and collection of annotation types,
#' the function queries the MGNify API and returns the results. This function
#' is convenient for retrieving highly structured (analysis vs counts) data on
#' certain instances. For example, BIOM files are downloaded automatically.
#' If you want just to retrieve raw data from the database, see \code{getData}.
#'
#' @param x A \code{MgnifyClient} object.
#'
#' @param accession A single character value or a vector of character values
#' specifying accession IDs to return results for.
#'
#' @param output A single character value specifying the format of an output.
#' Must be one of the following options: \code{"TreeSE"}, \code{"list"}, or
#' \code{"phyloseq"}. (By default: \code{output = "TreeSE"})
#'
#' @param get.taxa A boolean value specifying whether to retrieve taxonomy
#' data or character value specifying the type. If \code{bulk.dl=FALSE}, the
#' data is retrieved as BIOM files which are subsequently parsed.
#' (By default: \code{get.taxa = "SSU"})
#'
#' @param get.func A boolean value or a single character value or a vector
#' character values specifying functional analysis types to retrieve. If
#' \code{get.func = TRUE}, all available functional datatypes are retrieved,
#' and if \code{FALSE}, functional data is not retrieved. The current list of
#' available types is \code{"antismash-gene-clusters"}, \code{"go-slim"},
#' \code{"go-terms"}, and \code{"interpro-identifiers"}. Note that depending on
#' the particular analysis
#' type, pipeline version etc., not all functional results will be available.
#' Furthermore, taxonomy is also available via \code{get.func}, and loading
#' the data might be considerable faster if \code{bulk.dl = TRUE}. However,
#' phylogeny is available only via \code{get.taxa}.
#' (By default: \code{get.func = TRUE})
#'
#' @param ... optional arguments:
#' \itemize{
#'
#'   \item \strong{get.tree} A single boolean value specifying whether to
#'   include available phylogenetic trees in the \code{TreeSummarizedExperiment}
#'   object. Available when \code{get.taxa = TRUE}.
#'   (By default: \code{get.tree = TRUE})
#'
#'   \item \strong{as.df} A single boolean value enabled when
#'   \code{output = "list"}. The argument specifies whether return functional
#'   data as a named list (one entry per element in the output list) of
#'   data.frames, with each data.frame containing results for all requested
#'   accessions. If \code{FALSE}, the function returns a list of lists, each
#'   element consisting of results for a single accession. (By default:
#'   \code{as.df = TRUE})
#'
#'   \item \strong{bulk.dl} A single boolean value specifying should
#'   MGnifyR attempt to speed things up by downloading
#'   relevant studies TSV results and only extracting the required columns,
#'   rather than using the JSONAPI interface. When getting results where
#'   multiple accessions share the same study, this option may result in
#'   significantly faster processing. However, there appear to be (quite a few)
#'   cases in the database where the TSV result columns do NOT match the
#'   expected accession names. This will hopefully be fixed in the future,
#'   but for now \code{bulk.dl} defaults to TRUE. When it does work, it can
#'   be orders of magnitude more efficient.
#'   (By default: \code{buld.dl = TRUE})
#'
#' }
#'
#' @return
#' If only taxonomy data is retrieved, the result is returned in
#' \code{TreeSummarizedExperiment} object by default. The result can also be
#' returned as a \code{phyloseq} object or as a list of \code{data.frames}.
#' Note that \code{phyloseq} object can include only one phylogenetic tree
#' meaning that some taxa might be lost when data is subsetted based on tree.
#'
#' When functional data is retrieved in addition to taxonomy data, the result
#' is returned as a \code{MultiAssayExperiment} object. Other options are a list
#' containing \code{phyloseq} object and \code{data.frames} or just
#' \code{data.frames}.
#'
#' Functional data can be returned as a \code{MultiAssayExperiment} object or
#' as a list of \code{data.frames}.
#'
#' @examples
#' # Create a client object
#' mg <- MgnifyClient(useCache = FALSE)
#'
#' # Get OTU tables as TreeSE
#' accession_list <- c("MGYA00377505")
#' tse <- getResult(mg, accession_list, get.func=FALSE)
#'
#' \dontrun{
#' # Get functional data along with OTU tables as MAE
#' mae <- getResult(mg, accession_list, get.func=TRUE, get.taxa=TRUE)
#'
#' # Get same data as list
#' list <- getResult(
#'     mg, accession_list, get.func=TRUE, get.taxa=TRUE, output = "list",
#'     as.df = TRUE, use.cache = TRUE)
#' }
#'
#' @seealso
#' \code{\link[MGnifyR:getData]{getData}}
#'
#' @name getResult
NULL

#' @rdname getResult
#' @importFrom plyr llply
#' @importFrom dplyr bind_rows
#' @importFrom reshape2 dcast
#' @include AllClasses.R AllGenerics.R MgnifyClient.R utils.R
#' @export
setMethod("getResult", signature = c(x = "MgnifyClient"), function(
        x, accession, get.taxa = "SSU", get.func = TRUE, output = "TreeSE",
        ...){
    ############################### INPUT CHECK ################################
    if( !(.is_non_empty_character(accession)) ){
        stop(
            "'accession' must be a single character value or vector of ",
            "character values specifying the MGnify accession identifier.",
            call. = FALSE)
    }

    if( !(.is_a_bool(get.taxa) || is.character(get.taxa)) ){
        stop("'get.taxa' must be TRUE or FALSE or character value.",
            call. = FALSE)
    }
    if( !(.is_a_bool(get.func) || is.character(get.func)) ){
        stop("'get.func' must be TRUE or FALSE or character value.",
            call. = FALSE)
    }
    # If TRUE, fetch all the types
    all_types <- names(.analyses_results_type_parsers)
    supported_taxa <- all_types[ grepl("taxonomy", all_types) ]
    supported_func <- all_types[ !grepl("taxonomy", all_types) ]
    if( .is_a_bool(get.taxa) && get.taxa ){
        get.taxa <- supported_taxa
    }
    if( .is_a_bool(get.func) && get.func ){
        get.func <- supported_func
    }
    # Match the value with supported types
    if( is.character(get.taxa) ){
        get.taxa <- vapply(supported_taxa, function(x)
            grepl(paste0(get.taxa, collapse = "|"), x, ignore.case = TRUE),
            logical(1L))
        get.taxa <- get.taxa[ get.taxa ] |> names()
    }
    if( length(get.taxa) == 0L ){
        stop("'get.taxa' must be one of the following values: '",
            paste0(supported_taxa, collapse = "', '"), "'", call. = FALSE)
    }
    if( is.character(get.func) ){
        get.func <- vapply(supported_func, function(x)
            grepl(paste0(get.func, collapse = "|"), x, ignore.case = TRUE),
            logical(1L))
        get.func <- get.func[ get.func ] |> names()
    }

    if( length(get.func) == 0L ){
        stop("'get.func' must be one of the following values: '",
            paste0(supported_func, collapse = "', '"), "'", call. = FALSE)
    }
    if( !(length(output) == 1 && output %in% c("list", "phyloseq", "TreeSE")) ){
        stop(
            "'output' must be a 'TreeSE', 'list' or 'phyloseq'.", call. = FALSE)
    }
    ############################# INPUT CHECK END ##############################
    result <- .fetch_results(x, accession, get.func, get.taxa, output, ...)
    # Convert results to specified output type
    if( output != "list" ){
        result <- .convert_results_to_object(result, output)
    }
    return(result)
})

################################ HELP FUNCTIONS ################################

.fetch_results <- function(
        client, accession, get.func, get.taxa, output,
        use.cache = useCache(client), bulk.dl = FALSE, ...){
    if( !.is_a_bool(use.cache) ){
        stop("'use.cache' must be a single boolean value.", call. = FALSE)
    }
    if( !.is_a_bool(bulk.dl) ){
        stop("'bulk.dl' must be TRUE or FALSE.", call. = FALSE)
    }
    # If user specifies to fetch functional results from API instead of from
    # tsv files as bulk, the taxonomy data is available when fetching functional
    # data.
    if( bulk.dl ){
        get.func <- c(get.func, get.taxa)
        get.taxa <- FALSE
    }

    # Get functional data if user specified
    if( is.character(get.func) ){
        # If single value, create a vector from it
        if( length(get.func) == 1 ){
            get.func <- c(get.func)
        }
        func_res <- .mgnify_get_analyses_results(
            client = client, accession = accession, retrievelist = get.func,
            output = output, bulk.dl = bulk.dl, ...)
        # Put "taxonomy" first, then other taxonomy-related items, then the rest
        taxonomy_first <- intersect("taxonomy", names(func_res))
        taxonomy_related <- setdiff(
            grep("taxonomy|phylo-tax", names(func_res), value = TRUE),
            "taxonomy"
        )
        others <- setdiff(names(func_res), c(taxonomy_first, taxonomy_related))
        func_res <- func_res[c(taxonomy_first, taxonomy_related, others)]
        # In the old version, the "taxonomy" was called microbiota
        names(func_res)[ names(func_res) == "taxonomy" ] <- "microbiota"
    } else{
        func_res <- NULL
    }
    # Get microbial profiling data
    if( is.character(get.taxa) ){
        # The fetched BIOM files are parsed with mia::importBIOM, however,
        # mia does not import biomformat, it is only in its "suggests". This is
        # why we have to check that biomformat is available
        .require_package("biomformat")
        #
        taxa_res <- .mgnify_get_analyses_treese(
            client = client, accession = accession, get.taxa = get.taxa, ...)
    } else{
        taxa_res <- NULL
    }

    # The sample_data has been corrupted by doing the merge (names get messed
    # up and duplicated), so just regrab it with another lapply/rbind
    col_data <- lapply(accession, function(x){
        .mgnify_get_single_analysis_metadata(client, x, use.cache = use.cache)})

    result <- list(taxa = taxa_res, func = func_res, col_data = col_data)

    return(result)
}

# Convert results to TreeSE, MultiAssayExperiment or phyloseq.
#' @importFrom methods is
#' @importFrom dplyr bind_rows
.convert_results_to_object <- function(res, output){
    taxa_res <- res[["taxa"]]
    func_res <- res[["func"]]
    col_data <- res[["col_data"]]
    result <- NULL
    # If there are microbial profiling data, convert it to TreeSE or phyloseq
    if( !is.null(taxa_res) ){
        # Group based on analysis types
        tse_list <- .group_tse_list(taxa_res)
        # Merge data into TreeSE
        result <- lapply(tse_list, function(x){
            .merge_into_tse(x, col_data, output)
        })
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
                    .create_TreeSE_from_func_data,
                    tse = result[[1L]], col_data = col_data)
                # Get colData from microbial profiling data TreeSE,
                # if it is included
                args <- list()
                # If taxa data is included get all the data. Otherwise, get only
                # functional annotations.
                if( !is.null(result) ){
                    result <- c(result, func_res)
                } else {
                    result <- func_res
                }
            } else{
                # If user wants output as a phyloseq, give a list of one
                # phyloseq object and functional data
                result <- c(result, func_res)
                # If there are only one experiment, take it out from the list
                if( length(result) == 1 ){
                    result <- result[[1]]
                }
            }
        }
    }
    # If there are more than 1 experiments, create MAE
    if( output == "TreeSE" && length(result) > 1 ){
        col_data <- colData(result[[1L]])
        result <- ExperimentList(result)
        result <- MultiAssayExperiment(result)
        # If sample metadata is not empty
        if( !is.null(col_data) ){
            # Ensure that colData has all samples present in dataset
            all_samples <- unique(unlist(colnames(result)))
            col_data <- col_data[match(
                all_samples, rownames(col_data)), ]
            rownames(col_data) <- all_samples
            # And then add to mae
            colData(result) <- col_data
        }
    } else if( output == "TreeSE" ){
        # If there are only 1 experiment, give it as it is
        result <- result[[1]]
    }
    return(result)
}

# Convert a list of accession-type to type-accession
.group_tse_list <- function(tse_list){
    # get all unique taxonomic analysis names dynamically
    analysis_types <- unique(unlist(lapply(tse_list, names)))
    # Convert the nesting
    grouped <- lapply(analysis_types, function(type) {
        res <- lapply(names(tse_list), function(acc_name) {
            if (type %in% names(tse_list[[acc_name]])) {
                tse_list[[acc_name]][[type]]
            } else {
                NULL
            }
        })
        names(res) <- names(tse_list)
        res <- Filter(Negate(is.null), res)
        return(res)
    })
    names(grouped) <- analysis_types
    return(grouped)
}

.merge_into_tse <- function(tse_list, col_data, output){
    result <- NULL
    # If some results were not found, remove them
    ind <- !unlist(lapply(tse_list, is.null))
    # If there are samples left after subsetting
    if( any(ind) ){
        tse_list <- tse_list[ind]
        col_data <- col_data[ind]
        # Merge individual TreeSEs into one
        result <- mergeSEs(
            tse_list, assay.type = "counts", missing_values = 0)
        # Bind sample metadata to one table
        col_data <- do.call(bind_rows, col_data)
        col_data <- DataFrame(col_data)
        # Order the sample metadata
        col_data <- col_data[ colnames(result), , drop = FALSE]
        # Add sample metadata to the object
        colData(result) <- col_data
        # If user wants phyloseq, convert TreeSE
        if( output == "phyloseq" ){
            result <- makePhyloseqFromTreeSE(result)
        }
    } else{
        warning(
            "\nNo taxonomy data was found for the dataset.",
            call. = FALSE)
    }
    return(result)
}

# Create a TreeSE from single functional data data.frame
#' @importFrom S4Vectors SimpleList
#' @importFrom dplyr %>% mutate_all na_if
.create_TreeSE_from_func_data <- function(x, tse_microbiota, col_data){
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
        # Do not drop samples that are not found from microbiota data.
        # Order only when samples are shared.
        if( !is.null(tse_microbiota) &&
            all(colnames(assay) %in% colnames(tse_microbiota)) &&
            all(colnames(tse_microbiota) %in% colnames(assay)) ){
            assay <-  assay[ , match(
                colnames(tse_microbiota), colnames(assay)), drop = FALSE]
        }
        # Create a TreeSE
        assay <- as.matrix(assay)
        assays <- SimpleList(counts = assay)
        row_data <- DataFrame(row_data)
        # Get arguments for TreeSE
        args <- list(assays = assays, rowData = row_data)
        # Get sample metadata
        # Bind sample metadata to one table
        col_data <- do.call(bind_rows, col_data)
        col_data <- DataFrame(col_data)
        # Order coldata.
        col_data <- col_data[match(colnames(assay), rownames(col_data)), ]
        # Add colnames to ensure that all rows have name (if some samples
        # were missing from the col_data)
        rownames(col_data) <- colnames(assay)
        args$colData <- col_data
        # Create TreeSE
        tse <- do.call(TreeSummarizedExperiment, args)
        # If there are taxonomy ranks, the data includes taxonomy. Replace
        # names with tidier format.
        if( length(taxonomyRanks(tse)) > 1L ){
            rownames(tse) <- getTaxonomyLabels(tse)
        }
    } else{
        tse <- NULL
    }
    return(tse)
}

# Helper function for importing microbial profiling data.
.mgnify_get_analyses_treese <- function(
        client, accession, use.cache = useCache(client),
        show.messages = verbose(client), taxa.su = get.taxa,
        get.taxa = "SSU", ...){
    ############################### INPUT CHECK ################################
    if( !(.is_a_bool(taxa.su) || is.character(taxa.su)) ){
        stop(
            "'taxa.su' must be TRUE or FALSE.",
            call. = FALSE)
    }
    if( !.is_a_bool(use.cache) ){
        stop(
            "'use.cache' must be a single boolean value.", call. = FALSE)
    }
    if( !.is_a_bool(show.messages) ){
        stop(
            "'show.messages' must be a single boolean value.", call. = FALSE)
    }
    show.messages <- ifelse(show.messages, "text", "none")
    ############################# INPUT CHECK END ##############################
    # Get TreeSE objects
    names(accession) <- accession
    tse_list <- llply(accession, function(x) {
            .mgnify_get_single_analysis_treese(
                client, x, use.cache = use.cache, taxa.su = taxa.su, ...)
    }, .progress = show.messages)
    return(tse_list)
}

################################ HELP FUNCTIONS ################################

# Get a single biom file and convert it to TreeSummarizedExperiment format
#
#' @importFrom mia importBIOM
#' @importFrom mia checkTaxonomy
#' @importFrom urltools parameters parameters<-
#' @importFrom httr GET
#' @importFrom httr write_disk
#' @importFrom ape read.tree
#' @importFrom TreeSummarizedExperiment rowTree
#' @importFrom SummarizedExperiment rowData rowData<-
.mgnify_get_single_analysis_treese <- function(
        client = NULL, accession, use.cache = useCache(client),
        downloadDIR = NULL, taxa.su = "SSU", get.tree = TRUE,
        show.warnings = showWarnings(client), cache.dir = cacheDir(client),
        clear.cache = clearCache(client), ...){
    ############################### INPUT CHECK ################################
    if( !.is_a_bool(get.tree) ){
        stop("'get.tree' must be TRUE or FALSE.",
            call. = FALSE)
    }
    if( !.is_a_bool(use.cache) ){
        stop(
            "'use.cache' must be a single boolean value.", call. = FALSE)
    }
    if( !.is_a_bool(show.warnings) ){
        stop(
            "'show.warnings' must be a single boolean value.", call. = FALSE)
    }
    if( !.is_non_empty_string(cache.dir) ){
        stop(
            "'cache.dir' must be a string.", call. = FALSE)
    }
    if( !.is_a_bool(clear.cache) ){
        stop(
            "'clear.cache' must be a single boolean value.", call. = FALSE)
    }
    ############################# INPUT CHECK END ##############################
    # Get the metadata related to analysis
    metadata_df <- .mgnify_get_single_analysis_metadata(
        client, accession, use.cache = use.cache, ...)
    analysis_data <- .mgnify_retrieve_json(
        client, paste("analyses", accession, sep = "/"), use.cache = use.cache,
        ...)
    # Get the metadata related to available files etc
    download_url <- analysis_data[[1]]$relationships$downloads$links$related
    analysis_downloads <- .mgnify_retrieve_json(
        client, complete_url = download_url, use.cache = use.cache, ...)

    # Depending on the pipeline version, there may be more than one OTU table
    # available (LSU/SSU), so try and get the one specified in taxa.su -
    # otherwise spit out a warning and grab the generic (older pipelines)
    formats <- lapply(
        analysis_downloads, function(x){x$attributes$`file-format`$name})
    formats <- unlist(formats)
    available_biom_files <- analysis_downloads[grepl('JSON Biom', formats)]
    # Check if any biom files was found
    if( is.null(available_biom_files) || length(available_biom_files) == 0 ){
        warning(
            "\nNo BIOM data found for accession '", accession, "'.",
            call. = FALSE)
        return(NULL)
    }
    group_type <- lapply(
        available_biom_files, function(x){ x$attributes$`group-type` } )
    group_type <- unlist(group_type)
    names(available_biom_files) <- group_type
    taxa_to_fetch <- available_biom_files
    if( !.is_a_bool(taxa.su) ){
        pattern <- gsub("taxonomy-", "", paste0(taxa.su, collapse = "|"))
        biom_position <- grepl(pattern, group_type, ignore.case = TRUE)
        taxa_to_fetch <- available_biom_files[biom_position]
    }
    if( length(taxa_to_fetch) == 0L ){
        warning(
            "\nUnable to locate requested taxonomy type '", taxa.su, "'. ",
            "This is likely due to the current analysis having been ",
            "performed on an older version of the MGnify pipeline. ",
            "Use 'taxa.su' to specify the type. ",
            "The available BIOM files include '",
            paste0(group_type, collapse = "', '"), "'. Retrieving all of them.",
            call. = FALSE)
        taxa_to_fetch <- available_biom_files
    }
    biom_url <- lapply(
        taxa_to_fetch, function(x) x[["links"]][["self"]]) |> unlist()

    # Can specify a separate dir for saving biom files, otherwise they end up
    # in the cacheDir(client) folder, under "bioms"
    if( is.null(downloadDIR) ){
        downloadDIR <- file.path(cache.dir, "biom_files")
        dir.create(
            downloadDIR, recursive = TRUE, showWarnings = show.warnings)
    }

    # Fetch BIOM files
    tse_list <- lapply(biom_url, function(url){
        .fetch_single_biom_with_tree(
            biom_url = url, downloadDIR = downloadDIR, use.cache = use.cache,
            clear.cache = clear.cache, get.tree = get.tree,
            analysis_downloads = analysis_downloads, accession = accession)
    })
    return(tse_list)
}

# This function fetches a single BIOM file based on provided BIOM url and
# generates TreeSE from it.
.fetch_single_biom_with_tree <- function(
        biom_url, downloadDIR, use.cache, clear.cache, get.tree,
        analysis_downloads, accession){
    # Clear out any ?params after the main location - don't need them for this
    parameters(biom_url) <- NULL
    # Get the file name and path to it in local machine
    fname <- utils::tail(strsplit(biom_url, "/")[[1]], n = 1)
    biom_path <- file.path(downloadDIR, fname)

    ## Quick check to see if we should clear the disk cache  for this specific
    # call  - used for debugging and when MGnify breaks
    if( use.cache && clear.cache && file.exists(biom_path) ){
        message("clearCache is TRUE: deleting ", biom_path)
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
                "\n", biom_url, ": ", content(res, ...)$errors[[1]]$detail,
                " A biom listed in 'accession' is missing from the ",
                "output.", call. = FALSE)
            # Remove the downloaded file, it includes only info on errors
            unlink(biom_path)
            return(NULL)
        }
    }

    # Load in the TreeSummarizedExperiment object
    tse <- importBIOM(
        biom_path, rank.from.prefix = TRUE, set.ranks = TRUE, verbose = FALSE,
        prefix.rm = TRUE, artifact.rm = TRUE)
    # If the file was not in store already but fetched from database, and cache
    # storing is disabled
    if( fetched_from_url && !use.cache ){
        unlink(biom_path)
    }
    # TreeSE has sample ID as its colnames. Rename so that it is the accession
    # ID.
    colData(tse)[["biom_sample_id"]] <- colnames(tse)
    colnames(tse) <- accession

    # If user wants also phylogenetic tree
    if(get.tree){
        # Is there a tree?
        data_types <- lapply(analysis_downloads, function(x){
            x[["attributes"]][["description"]][["label"]]})
        data_types <- unlist(data_types)
        tvec <- grepl('Phylogenetic tree', data_types)
        if( any(tvec) ){
            # Get the url address of tree
            tree_url <- analysis_downloads[tvec][[1]]$links$self
            # Clear out any ?params after the main location - don't need them
            # for this
            parameters(tree_url) <- NULL
            # Get the file name for the tree
            fname <- utils::tail(strsplit(tree_url, '/')[[1]], n=1)
            # Get the path for the tree
            tree_path <- file.path(downloadDIR, fname)

            # Quick check to see if we should clear the disk cache
            #  for this specific call  - used for debugging and when MGnify
            # breaks
            if(use.cache && clear.cache && file.exists(tree_path) ){
                message("clearCache is TRUE: deleting ", tree_path)
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
                        "\n", tree_url, ": ",
                        content(res, ...)$errors[[1]]$detail,
                        " A phylogenetic tree listed in 'accession' is ",
                        "missing from the output.", call. = FALSE)
                }
            }
            # Add the tree to TreeSE object
            if( file.exists(tree_path) ){
                row_tree <- read.tree(tree_path)
                rowTree(tse) <- row_tree
                # If the file was not in store already but fetched from
                # database and cache storing is disabled
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
        client, accession, retrievelist, output, use.cache = useCache(client),
        show.messages = verbose(client), as.df = TRUE, bulk.dl = TRUE, ...){
    ############################### INPUT CHECK ################################
    # If output is TreeSE or object, get results as data.frames
    as.df <- ifelse(output != "list", TRUE, as.df)
    if( !.is_a_bool(as.df) ){
        stop("'as.df' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(bulk.dl) ){
        stop("'bulk.dl' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(use.cache) ){
        stop(
            "'use.cache' must be a single boolean value.", call. = FALSE)
    }
    if( !.is_a_bool(show.messages) ){
        stop(
            "'show.messages' must be a single boolean value.", call. = FALSE)
    }
    show.messages <- ifelse(show.messages, "text", "none")
    ############################# INPUT CHECK END ##############################
    # Get functional results
    all_results <- llply(accession, function(x){
        .mgnify_get_single_analysis_results(
            client, x, use.cache = use.cache, retrievelist = retrievelist,
            bulk.files = bulk.dl, ...)
    }, .progress = show.messages)
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
        client, accession, retrievelist=c(),
        use.cache = useCache(client), clear.cache = clearCache(client),
        max.hits = NULL, bulk.files = FALSE,
        show.warnings = showWarnings(client), ...){
    # Input check
    if( !.is_a_bool(clear.cache) ){
        stop(
            "'clear.cache' must be a single boolean value.", call. = FALSE)
    }
    if( !.is_a_bool(show.warnings) ){
        stop(
            "'show.warnings' must be a single boolean value.", call. = FALSE)
    }
    #
    # Get the metadata describing the samples
    metadata_df <- .mgnify_get_single_analysis_metadata(
        client, accession, use.cache = use.cache, max.hits = max.hits)

    # Should we try and grab the study's full TSV download rather than parse
    # through the JSON API? Doing so has the potential to use a LOT more disk
    # space, along with potentially increased data download. It should
    # be faster though, except in pathological cases (e.g. only 1 sample per
    # 1000 sample study required). As with everything else, we make use of
    # local caching to speed things along.
    if( bulk.files ){
        # Create a directory for donwloading the file
        downloadDIR <- file.path(cacheDir(client), "tsv")
        if(!dir.exists(downloadDIR)){
            dir.create(
                downloadDIR, recursive = TRUE, showWarnings = show.warnings)
        }
        # Get what files are available in database
        available_downloads_json <- .mgnify_retrieve_json(
            client,
            path = paste(
                "studies", metadata_df$study_accession, "downloads", sep = "/"),
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

            # Check the pipeline versions match and load the file if the
            # datatype is specified to be loaded.
            if(length(cur_type) > 0 &&
                cur_pipeversion == metadata_df$`analysis_pipeline-version` &&
                any(cur_type %in% retrievelist)){
                # If there are "taxonomic assignments ssu", it can match with
                # 2 --> take the first one
                cur_type <- cur_type[[1]]
                # Get the data
                temp <- .get_bulk_files(
                    cur_type, client = client, r = r,
                    metadata_df = metadata_df, downloadDIR = downloadDIR,
                    use.cache = use.cache, ...)
            } else{
                temp <- NULL
            }
            return(list(type=cur_type, data=temp))
        })
        # Add data types to data as names (taxonomy might have 2 data types if
        # NULL because these both 2 data types are tried to fetch)
        cur_type <- unlist(lapply(parsed_results, function(x)
            ifelse( length(x$type) > 1, x$type[[1]], x$type)))
        parsed_results <- lapply(parsed_results, function(x) x$data)
        names(parsed_results) <- cur_type
    }else{
        # If user do not want to fetch bulk files
        # Now (re)load the analysis data:
        analysis_data <- .mgnify_retrieve_json(
            client, paste("analyses", accession, sep = "/"),
            use.cache = use.cache, max.hits = max.hits, ...)
        # For now try and grab all data types that user has specified -
        # just return the list - don't do any processing...
        all_results <- lapply(
            names(.analyses_results_type_parsers), function(r) {
                tmp <- NULL
                if(r %in% retrievelist){
                    tmp <- .mgnify_retrieve_json(
                        client,
                        complete_url = analysis_data[[1]]$relationships[[
                            r]]$links$related,
                        use.cache = use.cache,
                        max.hits = max.hits, ...)
                }
                return(tmp)
        })
        # Add data types as names
        names(all_results) <- names(.analyses_results_type_parsers)
        # Get the data in correct format
        parsed_results <- lapply(names(all_results), function(x){
            # Get the specific type of data
            all_json <- all_results[[x]]
            # If specific type of data can be found
            res_df <- NULL
            if( !is.null(all_json) && length(all_json) > 0 ){
                # Get the specific type of data from all accessions
                all_json <- lapply(
                    all_json, .analyses_results_type_parsers[[x]])
                # Combine
                res_df <- do.call(bind_rows, all_json)
                rownames(res_df) <- res_df$index_id
            }
            return(res_df)
        })
        names(parsed_results) <- names(all_results)
    }
    # Return a list of data types, e.g., taxonomy, GO...
    return(parsed_results)
}

# Helper function for getting bulk file
.get_bulk_files <- function(
        cur_type, client, r, metadata_df, downloadDIR,
        use.cache, clear.cache = clearCache(client),
        use.mem.cache = FALSE, mem.cache.name = "mgnify_memory_cache", ...){
    # Input check
    if( !.is_a_bool(clear.cache) ){
        stop(
            "'clear.cache' must be a single boolean value.", call. = FALSE)
    }
    if( !.is_a_bool(use.mem.cache) ){
        stop(
            "'use.mem.cache' must be a single boolean value.", call. = FALSE)
    }
    # Check mem.cache.name
    if( !( length(mem.cache.name) == 1 && is.character(mem.cache.name) ) ){
        stop("'mem.cache.name' must be a single character value.",
            call. = FALSE)
    }
    if( exists(mem.cache.name) && use.mem.cache ){
        mgnify_memory_cache <- get(mem.cache.name)
        if( !is.list(mgnify_memory_cache) ){
            stop("'mem.cache.name' must specify a list object.", call. = FALSE)
        }
    } else if(use.mem.cache){
        mgnify_memory_cache <- list()
    }
    #
    # Get the url
    data_url <- r$links$self
    # Clear off extraneous gubbins
    parameters(data_url) <- NULL
    #build the cache filename
    fname <- utils::tail(strsplit(data_url, '/')[[1]], n=1)

    # At this point we might have already got the data we want
    # loaded. Check the memory cache object
    if( (use.mem.cache &&
        cur_type %in% names(mgnify_memory_cache) &&
        mgnify_memory_cache[cur_type]["fname"] == fname) ){
        tmp_df <- mgnify_memory_cache[cur_type][["data"]]
    }else{
        # Nope - gonna have to load it up from disk or grab
        # it from t'interweb
        data_path <- file.path(downloadDIR, fname)
        # Clear cache if specified
        if( use.cache && clear.cache && file.exists(data_path) ){
            message("clearCache is TRUE: deleting ", data_path)
            unlink(data_path)
        }
        # Download the file from the database to specific file path
        if( !file.exists(data_path) ){
            res <- GET(data_url, write_disk(data_path, overwrite = TRUE))
            # If the file was not successfully downloaded
            if( res$status_code != 200 ){
                warning(
                    "\n", data_url, ": ", content(res, ...)$errors[[1]]$detail,
                    " Could not load the file from database. The data ",
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
    if( use.mem.cache ){
        mgnify_memory_cache[[cur_type]] <- list(
            data=tmp_df, fname=fname)
        # Assign to global variable so it is visible for user
        assign(mem.cache.name, mgnify_memory_cache, envir = .GlobalEnv)
    }

    # Because there seem to be "mismatches" between the
    # JSON and downloadable files, and maybe some issues
    # with missing downloads, we have to check if we
    # actually got a valid file:
    if( ncol(tmp_df) < 3 ){
        warning("\nInvalid download for ", accession, call. = FALSE)
        return(NULL)
    }

    # Need to figure out how many columns to keep -
    # Get those columns that do not include numeric values but some info about
    # samples etc... First column includes always sample ID.
    info_cols <- lapply(tmp_df, function(x){
        all(is.na(suppressWarnings(as.numeric(x))))})
    info_cols <- unlist(info_cols)
    info_cols <- c(sample_ID = TRUE, info_cols[2:length(info_cols)])
    info_cols <- which(info_cols)

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
        warning("\nFailed to data on ", accession, call. = FALSE)
        return(NULL)
    }

    # Get the correct sample and subset the data so that it includes only
    # the sample and info columns
    column_position <- match(accession, colnames(tmp_df))
    if( (is.na(column_position) || length(column_position) != 1) ){
        warning("\nFailed to find column ", accession, call. = FALSE)
        return(NULL)
    }
    # Ensure that the counts are numeric
    tmp_df[[column_position]] <- as.numeric(tmp_df[[column_position]])
    # Get columns to keep
    keeper_columns <- c(info_cols, column_position)
    tmp_df <- tmp_df[, keeper_columns, drop = FALSE]
    # Adjust colnames rownames and add column
    colnames(tmp_df)[1] <- "accession"
    colnames(tmp_df)[ ncol(tmp_df) ] <- "count"
    tmp_df$index_id <- tmp_df$accession
    rownames(tmp_df) <- tmp_df$accession
    return(tmp_df)
}

# Which parser do you use for which type of output?
# a list of parsers for each output type. If new output types come, add them to
# here (and below) so that they are fetched from the database.
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
