#' Deprecated and defunct functions
#'
#' @details
#' MGnify API v1 is frozen, and MGnifyR now uses MGnify API v2. The functions
#' of the earlier versions of MGnifyR are deprecated, and they will be removed
#' in a future release. Where possible, they use the new functions:
#' \itemize{
#'   \item \code{MgnifyClient}: a client is not needed anymore. Log in with
#'   \code{\link{mgnifyLogin}} to access private data, and set
#'   \code{options(MGnifyR.cache = "path")} to keep the cache between
#'   sessions. For compatibility, \code{MgnifyClient} does these, and returns
#'   \code{NULL}. The \code{x} argument of the other deprecated functions is
#'   ignored.
#'   \item \code{getResult}: use \code{\link{fetchAnalyses}}. Only taxonomic
#'   profiles are returned, as \code{TreeSummarizedExperiment}.
#'   \item \code{searchAnalysis} and \code{getMetadata}: use
#'   \code{\link{searchAnalyses}}.
#'   \item \code{doQuery} and \code{getData}: use \code{\link{searchStudies}},
#'   \code{\link{searchSamples}} or \code{\link{searchAnalyses}}. These query
#'   any endpoint of
#'   \href{https://www.ebi.ac.uk/metagenomics/api/v2/}{MGnify API v2}, but
#'   the endpoints and filters differ from API v1.
#'   \item \code{searchFile} and \code{getFile}: list and download the files
#'   of studies, samples or analyses.
#' }
#' The accessors of \code{MgnifyClient}, such as \code{useCache}, and the
#' functions with \code{mgnify_} prefix are defunct.
#'
#' @param x Ignored. Formerly a \code{MgnifyClient} object.
#'
#' @param username,password Webin credentials for \code{\link{mgnifyLogin}}.
#'
#' @param useCache,cacheDir If \code{useCache = TRUE}, \code{cacheDir} is set
#' as the cache directory.
#'
#' @param accession A character vector of accessions.
#'
#' @param get.taxa Taxonomy to fetch, e.g., \code{"SSU"}.
#'
#' @param get.func,output Not supported anymore.
#'
#' @param type The type of accessions in \code{searchAnalysis} and
#' \code{searchFile}, or the API endpoint in \code{doQuery} and
#' \code{getData}, e.g., \code{"studies"}.
#'
#' @param accession.type The type of accessions in \code{getData}.
#'
#' @param as.df Whether to return a \code{data.frame} instead of a list.
#'
#' @param max.hits The maximum number of results, or \code{NULL} for all.
#'
#' @param url The URL of the file.
#'
#' @param file \code{NULL} or a path to save the file to.
#'
#' @param read.func \code{NULL} or a function applied to the downloaded file.
#'
#' @param value Ignored.
#'
#' @param ... In \code{doQuery}, query parameters of the API. Otherwise
#' ignored.
#'
#' @return See the replacing functions. \code{doQuery}, \code{getData} and
#' \code{searchFile} return a \code{data.frame}, and \code{getFile} returns
#' the path of the file or the result of \code{read.func}.
#'
#' @examples
#' # Instead of getResult(), use fetchAnalyses()
#' tse <- fetchAnalyses("MGYA00383350")
#'
#' @name deprecated
NULL

#' @rdname deprecated
#' @export
MgnifyClient <- function(
        username = NULL, password = NULL, useCache = FALSE,
        cacheDir = tempdir(), ...){
    .Deprecated(msg = paste0(
        "'MgnifyClient' is deprecated. A client is not needed anymore. Use ",
        "mgnifyLogin() to access private data, and options(MGnifyR.cache) ",
        "to set the cache directory."))
    if( !is.null(username) && !is.null(password) ){
        mgnifyLogin(username, password)
    }
    if( isTRUE(useCache) ){
        options(MGnifyR.cache = cacheDir)
    }
    return(invisible(NULL))
}

#' @rdname deprecated
#' @export
getResult <- function(
        x, accession, get.taxa = "SSU", get.func = TRUE, output = "TreeSE",
        ...){
    .Deprecated("fetchAnalyses")
    if( is.character(get.func) || !identical(output, "TreeSE") ){
        warning(
            "Only taxonomic profiles are returned as ",
            "TreeSummarizedExperiment.", call. = FALSE)
    }
    # Old taxonomies were, e.g., "taxonomy-ssu" or "taxonomy" for SSU
    type <- "SSU"
    if( is.character(get.taxa) ){
        type <- sub("^taxonomy-", "", get.taxa[[1L]])
    }
    if( !tolower(type) %in% tolower(.V2_TYPES[["type"]]) ){
        type <- "SSU"
    }
    res <- fetchAnalyses(accession, type = type)
    return(res)
}

#' @rdname deprecated
#' @export
searchAnalysis <- function(x, type, accession, ...){
    .Deprecated("searchAnalyses")
    if( !(length(type) == 1L && type %in% c("samples", "studies")) ){
        stop("'type' must be 'samples' or 'studies'.", call. = FALSE)
    }
    if( type == "studies" ){
        res <- searchAnalyses(accession)
        res <- setNames(res[["analysis_accession"]], res[["study_accession"]])
        return(res)
    }
    # Samples can be given with any of their accessions, so get the main
    # accession and the studies of the samples
    accession <- unique(accession)
    samples <- paste0("samples/", accession) |> .v2_url() |> .v2_get_json()
    .message_not_found(samples, accession, "Samples")
    names(samples) <- accession
    samples <- samples[lengths(samples) > 0L]
    studies <- lapply(samples, function(s){
        study_ids <- lapply(s[["studies"]], `[[`, "accession")
        return(study_ids)
    })
    studies <- unlist(studies) |> unique()
    res <- searchAnalyses(studies)
    ids <- vapply(samples, `[[`, character(1), "accession")
    res <- res[res[["sample_accession"]] %in% ids, , drop = FALSE]
    # Name the analyses by the sample accessions that were given
    rows <- match(res[["sample_accession"]], ids)
    given_ids <- names(ids)[rows]
    res <- setNames(res[["analysis_accession"]], given_ids)
    return(res)
}

#' @rdname deprecated
#' @export
getMetadata <- function(x, accession, ...){
    .Deprecated("searchAnalyses")
    # The old function returned all the metadata, so all the sources of
    # metadata are added
    res <- unique(accession) |> .v2_analyses_by_id(
        samples = TRUE, biosamples = TRUE, studies = TRUE, analyses = TRUE)
    rownames(res) <- res[["analysis_accession"]]
    return(res)
}

#' @rdname deprecated
#' @export
doQuery <- function(
        x, type = "studies", accession = NULL, as.df = TRUE, max.hits = 200,
        ...){
    .Deprecated("searchStudies")
    paths <- paste0(type, "/", accession)
    res <- .v2_query(paths, list(...), max.hits, as.df)
    return(res)
}

#' @rdname deprecated
#' @export
getData <- function(
        x, type, accession.type = NULL, accession = NULL, as.df = TRUE, ...){
    .Deprecated("searchAnalyses")
    paths <- paste0(type, "/")
    if( !is.null(accession.type) ){
        paths <- paste0(accession.type, "/", accession, "/", type, "/")
    }
    res <- .v2_query(paths, list(), NULL, as.df)
    return(res)
}

#' @rdname deprecated
#' @export
searchFile <- function(
        x, accession, type = c("studies", "samples", "analyses", "assemblies",
        "genomes", "run"), ...){
    .Deprecated(msg = "'searchFile' is deprecated.")
    type <- match.arg(type)
    accession <- unique(accession)
    items <- paste0(type, "/", accession) |> .v2_url() |> .v2_get_json()
    .message_not_found(items, accession, "Accessions")
    items <- items[lengths(items) > 0L]
    res <- lapply(items, function(x){
        df <- .v2_items_to_df(x[["downloads"]])
        if( nrow(df) > 0L ){
            df <- cbind(accession = x[["accession"]], type = type, df)
        }
        return(df)
    })
    # Start from an empty table, so that the result is a table even without
    # any files
    res <- c(list(data.frame()), res)
    res <- do.call(rbind, res)
    return(res)
}

#' @rdname deprecated
#' @export
getFile <- function(x, url, file = NULL, read.func = NULL, ...){
    .Deprecated(msg = "'getFile' is deprecated.")
    # Extension of the file, e.g., ".tsv", or "" if there is none
    ext <- tools::file_ext(url)
    ext <- if( nzchar(ext) ) paste0(".", ext) else ""
    path <- .v2_download(url, ext = ext)
    if( is.na(path) ){
        stop("Failed to download '", url, "'.", call. = FALSE)
    }
    if( !is.null(file) ){
        file.copy(path, file, overwrite = TRUE)
        path <- file
    }
    if( !is.null(read.func) ){
        res <- read.func(path)
        return(res)
    }
    return(path)
}

# Query any endpoint of the API. List endpoints are fetched with pagination,
# and detail endpoints as they are.
.v2_query <- function(paths, query, max.hits, as.df){
    res <- .v2_url(paths, query) |> .v2_get_json()
    .message_not_found(res, paths, "Data")
    items <- lapply(seq_along(paths), function(i){
        # Detail endpoint, e.g., "studies/MGYS1", has a single item
        if( is.null(res[[i]][["items"]]) ){
            item <- res[i]
            return(item)
        }
        # List endpoint, e.g., "studies/", has many pages of items
        path_items <- .v2_get_all(paths[[i]], query, max.hits)[[1L]]
        return(path_items)
    })
    items <- unlist(items, recursive = FALSE)
    items <- items[lengths(items) > 0L]
    if( as.df ){
        items <- .v2_items_to_df(items)
    }
    return(items)
}

#' @rdname deprecated
#' @export
databaseUrl <- function(...) .defunct_client("databaseUrl")
#' @rdname deprecated
#' @export
authTok <- function(...) .defunct_client("authTok")
#' @rdname deprecated
#' @export
useCache <- function(...) .defunct_client("useCache")
#' @rdname deprecated
#' @export
cacheDir <- function(...) .defunct_client("cacheDir")
#' @rdname deprecated
#' @export
showWarnings <- function(...) .defunct_client("showWarnings")
#' @rdname deprecated
#' @export
clearCache <- function(...) .defunct_client("clearCache")
#' @rdname deprecated
#' @export
verbose <- function(...) .defunct_client("verbose")
#' @rdname deprecated
#' @export
`databaseUrl<-` <- function(x, value) .defunct_client("databaseUrl<-")
#' @rdname deprecated
#' @export
`authTok<-` <- function(x, value) .defunct_client("authTok<-")
#' @rdname deprecated
#' @export
`useCache<-` <- function(x, value) .defunct_client("useCache<-")
#' @rdname deprecated
#' @export
`cacheDir<-` <- function(x, value) .defunct_client("cacheDir<-")
#' @rdname deprecated
#' @export
`showWarnings<-` <- function(x, value) .defunct_client("showWarnings<-")
#' @rdname deprecated
#' @export
`clearCache<-` <- function(x, value) .defunct_client("clearCache<-")
#' @rdname deprecated
#' @export
`verbose<-` <- function(x, value) .defunct_client("verbose<-")

.defunct_client <- function(name){
    .Defunct(msg = paste0(
        "'", name, "' is defunct. A client is not needed anymore. Use ",
        "mgnifyLogin() to access private data, and options(MGnifyR.cache) ",
        "to set the cache directory."))
}

#' @rdname deprecated
#' @export
mgnify_client <- function(...) .Defunct("MgnifyClient")
#' @rdname deprecated
#' @export
mgnify_query <- function(...) .Defunct("searchStudies")
#' @rdname deprecated
#' @export
mgnify_analyses_from_samples <- function(...) .Defunct("searchAnalyses")
#' @rdname deprecated
#' @export
mgnify_analyses_from_studies <- function(...) .Defunct("searchAnalyses")
#' @rdname deprecated
#' @export
mgnify_get_download_urls <- function(...) .Defunct("searchFile")
#' @rdname deprecated
#' @export
mgnify_download <- function(...) .Defunct("getFile")
#' @rdname deprecated
#' @export
mgnify_get_analyses_results <- function(...) .Defunct("fetchAnalyses")
#' @rdname deprecated
#' @export
mgnify_get_analyses_phyloseq <- function(...) .Defunct("fetchAnalyses")
#' @rdname deprecated
#' @export
mgnify_get_analyses_metadata <- function(...) .Defunct("searchAnalyses")
#' @rdname deprecated
#' @export
mgnify_retrieve_json <- function(...) .Defunct("searchStudies")
