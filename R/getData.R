#' Versatile function to retrieve raw results
#'
#' @details
#' This function returns data from MGnify database. Compared to
#' \code{getResult}, this function allows more flexible framework for fetching
#' the data. However, there are drawbacks: for counts data, \code{getResult}
#' returns optimally structured data container which is easier for downstream
#' analysis. \code{getData} returns raw data from the database. However, if
#' you want to retrieve data on pipelines or publications, for instance,
#' \code{getResult} is not suitable for it, and \code{getData} can be utilized
#' instead.
#'
#' @param x A \code{MgnifyClient} object.
#'
#' @param type A single character value specifying the type of data retrieve.
#' Must be one of the following options: \code{studies}, \code{samples},
#' \code{runs}, \code{analyses}, \code{biomes}, \code{assemblies},
#' \code{super-studies}, \code{experiment-types}, \code{pipelines},
#' \code{pipeline-tools}, \code{publications}, \code{genomes},
#' \code{genome-search}, \code{genome-search/gather}, \code{genome-catalogues},
#' \code{genomeset}, \code{cogs}, \code{kegg-modules}, \code{kegg-classes},
#' \code{antismash-geneclusters}, \code{annotations/go-terms},
#' \code{annotations/interpro-identifiers}, \code{annotations/kegg-modules},
#' \code{annotations/pfam-entries}, \code{annotations/kegg-orthologs},
#' \code{annotations/genome-properties},
#' \code{annotations/antismash-gene-clusters}, \code{annotations/organisms}, or
#' \code{mydata}.
#'
#' @param accession A single character value or a vector of character values
#' specifying accession IDs to return results for.
#' (By default: \code{accession = NULL})
#'
#' @param accession.type A single character value specifying type of accession
#' IDs (\code{accession}). Must be specified when \code{accession} is specified.
#' (By default: \code{accession.type = NULL})
#'
#' @param as.df A single boolean value specifying whether to return the
#' results as a data.frame or leave as a nested list.
#' (By default: \code{as.df = TRUE})
#'
#' @param ... optional arguments fed to internal functions.
#'
#' @return
#' \code{data.frame} or \code{list}
#'
#' @examples
#' # Create a client object
#' mg <- MgnifyClient(useCache = FALSE)
#'
#' # Find kegg modules for certain analysis
#' df <- getData(
#'     mg, type = "kegg-modules",
#'     accession = "MGYA00642773", accession.type = "analyses")
#'
#' @seealso
#' \code{\link[MGnifyR:getResult]{getResult}}
#'
#' @name getData
NULL

#' @rdname getData
#' @include allClasses.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "getData", signature = c(x = "MgnifyClient"), function(
    x, type, accession.type = NULL, accession = NULL, as.df = TRUE, ...){
    ############################### INPUT CHECK ################################
    available_types <- c(
        "studies", "samples", "runs", "analyses", "biomes", "assemblies",
        "super-studies", "experiment-types", "pipelines", "pipeline-tools",
        "publications", "genomes", "genome-search", "genome-search/gather",
        "genome-catalogues", "genomeset", "cogs", "kegg-modules",
        "kegg-classes", "antismash-geneclusters", "annotations/go-terms",
        "annotations/interpro-identifiers", "annotations/kegg-modules",
        "annotations/pfam-entries", "annotations/kegg-orthologs",
        "annotations/genome-properties", "annotations/antismash-gene-clusters",
        "annotations/organisms", "mydata")
    if( !(.is_non_empty_string(type) && type %in% available_types) ){
        stop(
            "'type' must be a single character value specifying ",
            "the type of instance to query. The value must be one of the ",
            "following options: ",
            paste0("'", paste(available_types, collapse = "', '"), "'"),
            call. = FALSE)
    }
    if( !(.is_non_empty_character(accession) || is.null(accession)) ){
        stop(
            "'accession' must be a single character value or vector of ",
            "character values specifying the MGnify accession identifier.",
            call. = FALSE)
    }
    if( !(.is_non_empty_character(accession.type) || is.null(accession.type)) ){
        stop(
            "'accession.type' must be a single character value or vector of ",
            "character values specifying the type of MGnify accession ",
            "identifier.", call. = FALSE)
    }
    if(
        (is.null(accession) && !is.null(accession.type)) ||
        (is.null(accession.type) && !is.null(accession)) ){
        stop(
          "Both 'accession' and 'accession.type' must be specified or they ",
          "must be NULL.", call. = FALSE)
    }
    if( !.is_a_bool(as.df) ){
        stop(
            "'as.df' must be a single boolean value specifying whether",
            "to return list or data.frame.", call. = FALSE)
    }
    ############################# INPUT CHECK END ##############################
    # Retrieve results
    result <- .get_results_as_json_list(x, type, accession.type, accession, ...)
    # Convert to df
    if( as.df ){
        result <- .convert_json_list_to_df(result)
    } else if( length(result) == 1 ){
        result <- result[[1]]
    }
    return(result)
})

################################ HELP FUNCTIONS ################################

#' @importFrom plyr llply
.get_results_as_json_list <- function(mg, type, accession.type, accession, ...){
    # Create a path. If multiple accession IDs, path is vector of multiple
    # paths. Otherwise the path specifies only the type
    if( !is.null(accession.type) && !is.null(accession) ){
        path <- paste0(accession.type, "/", accession, "/", type)
        names(path) <- accession
    } else{
        path <- type
    }
    # Find results by loping through paths
    res <- llply(path, function(x){
        .mgnify_retrieve_json(mg, path = x, ...)
    })
    return(res)
}

#' @importFrom tidyjson spread_all
#' @importFrom dplyr bind_rows
.convert_json_list_to_df <- function(result){
    # Create data.frames from individual search results
    res <- lapply(result, function(x){
        if( !is.null(x) ){
            x <- as.data.frame(spread_all(x))
        }
        return(x)
    })
    # Merge individual data.frames to one
    res <- bind_rows(res)
    # Add names if there were accession IDs provided as input
    if( !is.null(names(result)) ){
        # Assign to "accession" column name if there is no column with that name
        # already
        col_name <- "accession"
        col_name <- c(colnames(res), col_name)
        col_name <- make.unique(col_name)[[ length(col_name) ]]
        # Add to result df
        nams <- rep( names(result), each = lengths(result))
        res[[ col_name ]] <- nams
    }
    return(res)
}
