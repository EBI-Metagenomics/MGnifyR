#' These functions will be deprecated. Please use other functions instead.
#'
#' @param url -
#'
#' @param username -
#'
#' @param password -
#'
#' @param usecache -
#'
#' @param cache_dir -
#'
#' @param warnings -
#'
#' @param use_memcache -
#'
#' @param client -
#'
#' @param qtype -
#'
#' @param accession -
#'
#' @param asDataFrame -
#'
#' @param maxhits -
#'
#' @param ... -
#'
#' @param accessions -
#'
#' @param accession_type -
#'
#' @param file -
#'
#' @param read_func -
#'
#' @param Debug -
#'
#' @param retrievelist -
#'
#' @param compact_results -
#'
#' @param bulk_dl -
#'
#' @param returnLists -
#'
#' @param tax_SU -
#'
#' @param get_tree -
#'
#' @param path -
#'
#' @param complete_url -
#'
#' @param qopts -
#'
#' @name deprecate
NULL

#' @rdname deprecate
#' @export
mgnify_client <- function(
        url = NULL, username = NULL, password = NULL, usecache = FALSE,
        cache_dir = NULL, warnings = FALSE, use_memcache = FALSE){
    .Deprecated("MgnifyClient")
    MgnifyClient(url = url,
                 username = username, password = password,
                 use.cache = usecache, cache.dir = cache_dir, warnings = warnings,
                 use.memcache = use_memcache)
}

#' @rdname deprecate
#' @export
mgnify_query <- function(
        client, qtype = "samples", accession = NULL, asDataFrame = TRUE,
        maxhits = 200, usecache = FALSE, ...){
    .Deprecated("doQuery")
    doQuery(
        x = client, type = qtype, accession = accession,
        as.df = asDataFrame, max.hits = maxhits, usecache = usecache, ...)
}

#' @rdname deprecate
#' @export
mgnify_analyses_from_samples <- function(
        client, accession, usecache = TRUE, ...){
    .Deprecated("getAnalysisAccessions")
    getAnalysisAccessions(
        x = client, type = "sample", accession = accession,
        use.cache = usecache, ...)
}

#' @rdname deprecate
#' @export
mgnify_analyses_from_studies <- function(
        client, accession, usecache = TRUE, ...){
    .Deprecated("getAnalysisAccessions")
    getAnalysisAccessions(
        x = client, type = "study", accession = accession,
        use.cache = usecache, ...)
}

#' @rdname deprecate
#' @export
mgnify_get_download_urls <- function(
        client, accessions, accession_type, usecache = TRUE, ...){
    .Deprecated("getFileUrl")
    getFileUrl(
        x = client, accession = accessions, type = accession_type,
        use.cache = usecache, ...)
}

#' @rdname deprecate
#' @export
mgnify_download <- function(
        client, url, file = NULL, read_func = NULL, usecache = TRUE,
        Debug = FALSE, ...){
    .Deprecated("getFile")
    getFile(
        x = client, url = url, file = file,
        read.func = read_func, use.cache = usecache, ...)
}

#' @rdname deprecate
#' @export
mgnify_get_analyses_results <- function(
        client=NULL, accessions, retrievelist = c(), compact_results = TRUE,
        usecache = TRUE, bulk_dl = FALSE, ...){
    .Deprecated("getResults")
    if( length(retrievelist) == 0 ){
        retrievelist <- FALSE
    }
    getResults(
        x = client, accession = accessions, get.taxa = FALSE,
        get.func = retrievelist, output = "list", usecache = TRUE,
        as.df = compact_results, ...)
}

#' @rdname deprecate
#' @export
mgnify_get_analyses_phyloseq <- function(
        client = NULL, accessions, usecache = TRUE, returnLists = FALSE,
        tax_SU = "SSU", get_tree = FALSE, ...){
    .Deprecated("getResults")
    output <- ifelse(returnLists, "list", "phyloseq")
    getResults(
        x = client, accession = accessions, get.taxa = TRUE, get.func = FALSE,
        output = output, use.cache = usecache, tax.su = tax_SU,
        get.tree = get_tree, ...
    )
}

#' @rdname deprecate
#' @export
mgnify_get_analyses_metadata <- function(
        client, accessions, usecache = TRUE, ...){
    .Deprecated("getMetadata")
    getMetadata(x = client, accession = accessions, usecache = usecache, ...)
}

#' @rdname deprecate
#' @export
mgnify_retrieve_json <- function(
        client, path = "biomes", complete_url = NULL, qopts = NULL,
        maxhits = 200, usecache = FALSE, Debug = FALSE){
    .Deprecated(msg = "'mgnify_retrieve_json' is deprecated.\n",
                "See other functions and use them instead.\n",
                "See help('Deprecated')")
    return(NULL)
}
