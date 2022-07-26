#' Get functional and/or taxonomic information for a list of accessions
#'
#' Given a set of analysis accessions and collection of annotation types, \code{mgnify_get_analyses_results} queries the MGNify API
#' and returns the results, by default merging the results into multi-accession data.frames
#'
#' @importFrom plyr llply
#' @importFrom dplyr bind_rows
#' @importFrom reshape2 dcast
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
#' @return Named list of \code{data.frames}, corresponding to the requested analysis types in \code{retrievelist}
#' @examples
#'
#'@export
mgnify_get_analyses_results <- function(client=NULL, accessions, retrievelist=c(), compact_results=T, usecache = T, bulk_dl = F){
    if(length(retrievelist) == 1 && retrievelist == "all"){
        retrievelist <- names(analyses_results_type_parsers)
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
                r <- lapply(results_as_lists, function(x){
                    df <- as.data.frame(x[[y]])
                    df
                })
                longform <- dplyr::bind_rows(r, .id = "analysis")
                cn <- colnames(longform)
                extras <- cn[!(cn %in% c("count","index_id", "analysis"))]
                final_df <- reshape2::dcast(longform, as.formula(paste(paste(extras,collapse = " + "), " ~ analysis")), value.var = "count", fun.aggregate = sum)
                final_df}, error=function(x) NULL)
        })
    }
    names(all_results) <- retrievelist
    all_results
}
