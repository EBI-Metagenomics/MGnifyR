#' Retrieve OTU tables for all specified accessions and build a\code{phyloseq} object
#'
#' \code{mgnify_get_analyses_phyloseq} retrieves all associated Study, Sample and Analysis metadata attributes,
#' along with all OTU tables (of a given taxonomic type), and merges them together to build a \code{phyloseq}
#' object with \code{otu_table}, \code{tax_table} and \code{sample_data} objects.
#'
#' @importFrom plyr llply
#' @importFrom phyloseq merge_phyloseq
#' @importFrom phyloseq sample_data
#' @importFrom dplyr bind_rows
#'
#' @param client \code{mgnify_client} instance
#' @param accessions Single value or list/vector of Analysis accessions to retrieve data for.
#' @param usecache Whether to use the disk based cache.
#' @param returnLists Flag to determine whether to merge a per-analysis \code{phyloseq} objects into one, or
#' return a list of single-sample objects. Since perculiarities in OTU format may lead to failed merges, returning
#' a list of single objects may aid debugging. In most cases though the desired behaviour is to return a single
#' object (returnLists = F)
#' @param tax_SU Which taxa subunit results should be selected? Currently, taxonomy assignments in the
#' MGnify pipelines rely on rRNA matches to existing databases (GreenGenes and SILVA), with later
#' pipelines checking both the SSU and LSU portions of the rRNA sequence. \code{tax_SU} allows the
#' selection of either the Small subunit (SSU) or Large subunit results in the final \code{phyloseq} object.
#' Older pipeline versions do not report results for both subunits,
#' and thus for some accessions this value will have no effect.
#' @param get_tree Flag to control whether to include available phylogenetic trees in the phyloseq object. At present this option is of limited use as
#' the trees available are specific to each accession, holding only a subset of the leaves and branches of the full canonical tree used for the pipeline. This means that
#' \code{phyloseq} is unable to merge the trees, and therefore fails to build a final combined object. Setting \code{returnLists} to TRUE allows the results to be returned successfully,
#' albeit not in a single object.
#' @return Combined phyloseq object with \code{otu_table}, \code{sample_data} and \code{tax_table} entries for all accessions.
#' @examples
#'
#'
#' @export
mgnify_get_analyses_phyloseq <- function(client = NULL, accessions, usecache=T,
                                                                                 returnLists=F, tax_SU = "SSU",
                                                                                 get_tree=FALSE){
    #Some biom files don't import - so need a try/catch
    ps_list <- llply(accessions, function(x) {
        tryCatch(
                mgnify_get_single_analysis_phyloseq(client, x, usecache = usecache, tax_SU = tax_SU, get_tree = get_tree), error=function(x) NULL)
        }, .progress = "text")

    #The sample_data has been corrupted by doing the merge (names get messed up and duplicated), so just regrab it with another lapply/rbind
    samp_dat <- lapply(accessions, function(x) mgnify_get_single_analysis_metadata(client, x, usecache = usecache ))
    if (returnLists){
        list(phyloseq_objects=ps_list, sample_metadata = samp_dat)
    }else{

        #first of all, check to see that if we wanted them, we got trees for ALL the phyloseq objects.
        #If trees are present in any of the phyloseq objects during merging, then any OTUs not in a tree
        #(e.g. if any phyloseq objects do NOT contain a tree) will not be included in the merged output.

        if(get_tree){
            if (any(is.na(lapply(ps_list, function(x) x@phy_tree)))){
            warning("Phylogenetic tree retrieval was requested but some of the analyses do not include phylogenetic trees. Results should be used with caution.")
            }
        }

        #This is too slow for large datasets
        #full_ps <- do.call(merge_phyloseq, ps_list)
        #so:
        #a divide-and-conquer approach to merge_phyloseq seems to work best, hence the following
        #code which splits the full list into sublists and merges them seperately, then repeats until all are joined.
        curlist <- ps_list
        while(length(curlist) > 1){
            #Lists of length 10 seem to work well
            sublist <- split(curlist, seq_along(curlist) %/% 10)
            curlist <- lapply(sublist, function(x){
                do.call(merge_phyloseq,x)
            })
        }
        #By this point curlist isn't a list, it's a phyloseq object...
        full_ps <- curlist[[1]]
        sample_metadata_df <- do.call(bind_rows, samp_dat)
        rownames(sample_metadata_df) <- sample_metadata_df$analysis_accession
        sample_data(full_ps) <- sample_metadata_df
        full_ps
    }
}
