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
