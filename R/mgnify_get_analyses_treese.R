#' Retrieve OTU tables for all specified accessions and build a\code{TreeSummarizedExperiment} object
#'
#' \code{mgnify_get_analyses_treese} retrieves all associated Study, Sample and Analysis metadata attributes,
#' along with all OTU tables (of a given taxonomic type), and merges them together to build a \code{TreeSummarizedExperiment}
#' object with \code{otu_table}, \code{tax_table} and \code{sample_data} objects.
#'
#' @importFrom plyr llply
#' @importFrom mia mergeSEs
#'
#' @param client \code{mgnify_client} instance
#' @param accessions Single value or list/vector of Analysis accessions to retrieve data for.
#' @param usecache Whether to use the disk based cache.
#' @param returnLists Flag to determine whether to merge a per-analysis \code{TreeSummarizedExperiment} objects into one, or
#' return a list of single-sample objects. Since perculiarities in OTU format may lead to failed merges, returning
#' a list of single objects may aid debugging. In most cases though the desired behaviour is to return a single
#' object (returnLists = F)
#' @param tax_SU Which taxa subunit results should be selected? Currently, taxonomy assignments in the
#' MGnify pipelines rely on rRNA matches to existing databases (GreenGenes and SILVA), with later
#' pipelines checking both the SSU and LSU portions of the rRNA sequence. \code{tax_SU} allows the
#' selection of either the Small subunit (SSU) or Large subunit results in the final \code{TreeSummarizedExperiment} object.
#' Older pipeline versions do not report results for both subunits,
#' and thus for some accessions this value will have no effect.
#' @param get_tree Flag to control whether to include available phylogenetic trees in the TreeSummarizedExperiment object. At present this option is of limited use as
#' the trees available are specific to each accession, holding only a subset of the leaves and branches of the full canonical tree used for the pipeline. This means that
#' \code{TreeSummarizedExperiment} is unable to merge the trees, and therefore fails to build a final combined object. Setting \code{returnLists} to TRUE allows the results to be returned successfully,
#' albeit not in a single object.
#' @return Combined TreeSummarizedExperiment object with \code{row_data}, \code{col_data} and \code{assays} entries for all accessions.
#' @examples
#'
#'
#' @export
mgnify_get_analyses_treese <- function(client = NULL, accessions, usecache=T, returnLists=F, tax_SU = "SSU", get_tree=FALSE){
    #Some biom files don't import - so need a try/catch
    tse_list <- plyr::llply(accessions, function(x) {
        tryCatch(
            mgnify_get_single_analysis_treese(client, x, usecache = usecache, tax_SU = tax_SU, get_tree = get_tree),
            warning=function(x){ message("a biom listed in \"accessions\" is missing from the retrieved tse_list")},
            error=function(x) NULL
        )
    }, .progress = "text")

    #The sample_data has been corrupted by doing the merge (names get messed up and duplicated), so just regrab it with another lapply/rbind
    samp_dat <- lapply(accessions, function(x) mgnify_get_single_analysis_metadata(client, x, usecache = usecache ))
    if (returnLists){
        list(tse_objects=tse_list, sample_metadata = samp_dat)
    }else{
        if(get_tree){
            if (any(is.na(lapply(tse_list, function(x) x@rowTree)))){
                warning("Phylogenetic tree retrieval was requested but some of the analyses do not include phylogenetic trees. Results should be used with caution.")
            }
        }
        
        #Merging every TSE object in the list into a single one TSE object
        full_tse <- mia::mergeSEs(tse_list, assay_name = "counts", missing_values = 0)
        full_tse
    }
}
