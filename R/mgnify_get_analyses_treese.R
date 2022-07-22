mgnify_get_analyses_treese <- function(client = NULL, accessions, usecache=T,
                                    returnLists=F, tax_SU = "SSU",
                                    get_tree=FALSE){
  #Some biom files don't import - so need a try/catch
  tse_list <- plyr::llply(accessions, function(x) {
    tryCatch(
      mgnify_get_single_analysis_treese(client, x, usecache = usecache, tax_SU = tax_SU, get_tree = get_tree), error=function(x) NULL)
  }, .progress = "text")

  #The sample_data has been corrupted by doing the merge (names get messed up and duplicated), so just regrab it with another lapply/rbind
  samp_dat <- lapply(accessions, function(x) mgnify_get_single_analysis_metadata(client, x, usecache = usecache ))
  if (returnLists){
    list(tse_objects=tse_list, sample_metadata = samp_dat)
  }else{

    if(get_tree){
      if (any(is.na(lapply(ps_list, function(x) x@rowTree)))){
        warning("Phylogenetic tree retrieval was requested but some of the analyses do not include phylogenetic trees. Results should be used with caution.")
      }
    }

    full_tse <-mia::mergeSEs(tse_list, assay_name = "abundance", missing_values = 0)
    full_tse

  }
}
