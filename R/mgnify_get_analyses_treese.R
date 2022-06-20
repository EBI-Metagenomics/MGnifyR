mgnify_get_single_analysis_treese <- function(client=NULL, accession, usecache=T, downloadDIR=NULL, tax_SU="SSU", get_tree=FALSE){

  metadata_df <- mgnify_get_single_analysis_metadata(client, accession, usecache=usecache)
  analysis_data <-  mgnify_retrieve_json(client, paste("analyses",accession,sep="/"), usecache = usecache)
  download_url <- analysis_data[[1]]$relationships$downloads$links$related
  analysis_downloads <- mgnify_retrieve_json(client, complete_url = download_url,usecache = usecache)

  #Depending on the pipeline version, there may be more than one OTU table available (LSU/SSU), so try and get the
  #one specified in tax_SU - otherwise spit out a warning and grab the generic (older pipelines)
  available_biom_files <- analysis_downloads[grepl('JSON Biom', sapply(analysis_downloads, function(x) x$attributes$`file-format`$name))]
  biom_position <- grepl(tax_SU, sapply(available_biom_files, function(x) x$attributes$`group-type`))
  if(sum(biom_position) == 0){
    if(client@warnings){
    warning("Unable to locate requested taxonomy type ",tax_SU,". This is likely due to the current analysis having been performed on an older version of the MGnify pipeline.
             The available BIOM file will be used instead.")
    }
    biom_url <- available_biom_files[[1]]$links$self
  }else{
    biom_url <- available_biom_files[biom_position][[1]]$links$self
  }

  #Can specify a seperate dir for saving biom files, otherwise they end up in the client@cachdir folder, under "bioms"
  if (is.null(downloadDIR)){
    downloadDIR=paste(client@cache_dir,"biom_files",sep="/")
    dir.create(downloadDIR, recursive = T, showWarnings = client@warnings)
  }
  #Clear out any ?params after the main location - don't need them for this
  urltools::parameters(biom_url) <- NULL

  fname=tail(strsplit(biom_url, '/')[[1]], n=1)
  biom_path = paste(downloadDIR, fname, sep="/")

  ## Quick check to see if we should clear the disk cache ~for this specific call~ - used for debugging
  # and when MGnify breaks
  if(usecache & client@clear_cache){
    print(paste("clear_cache is TRUE: deleting ",biom_path, sep=""))
    tryCatch(unlink(biom_path), error=warning)
  }

  if (! file.exists(biom_path)){#} | !use_downloads ){
    httr::GET(biom_url, httr::write_disk(biom_path, overwrite = T))
  }

  #Load in a SummarizedExperiment object
  tse <- mia::loadFromBiom(biom_path)
  #Need to check if the taxonomy was parsed correctly - depending on the pipeline it may need a bit of help:
  mia::checkTaxonomy(tse)

  #Turn the SummarizedExperiment object into a TreeSummarizedExperiment object
  tse <- as(se, "TreeSummarizedExperiment")

  #The biom files have a single column of unknown name - maybe it's the original sample name?.
  # It's rewritten as sample_run_analysis accession, with
  # the original value stored in the sample_data (just in case it's needed later)
  orig_samp_name <- colnames(tse)
  metadata_df[1,"orig_samp_name"] <- orig_samp_name
  colData(tse) <- metadata_df

  if(get_tree){
    #is there a tree?
    tvec = grepl('Phylogenetic tree', sapply(analysis_downloads, function(x) x$attributes$`description`$label))
    if(any(tvec)){
      tree_url = analysis_downloads[tvec][[1]]$links$self
      #Clear out any ?params after the main location - don't need them for this
      urltools::parameters(tree_url) <- NULL

      fname=tail(strsplit(tree_url, '/')[[1]], n=1)
      tree_path = paste(downloadDIR, fname, sep="/")

      ## Quick check to see if we should clear the disk cache ~for this specific call~ - used for debugging
      # and when MGnify breaks
      if(usecache & client@clear_cache){
        print(paste("clear_cache is TRUE: deleting ",tree_path, sep=""))
        tryCatch(unlink(tree_path), error=warning)
      }

      if (! file.exists(tree_path)){#} | !use_downloads ){
        httr::GET(tree_url, httr::write_disk(tree_path, overwrite = T ))
      }
    }
    row_tree = ape::read.tree(tree_path)
    rowTree(tse) <- row_tree
  }
  tse
}


mgnify_get_analyses_treese <- function(client = NULL, accessions, usecache=T,
                                    returnLists=F, tax_SU = "SSU",
                                    get_tree=FALSE){
  #Some biom files don't import - so need a try/catch
  ps_list <- plyr::llply(accessions, function(x) {
    tryCatch(
      mgnify_get_single_analysis_treese(client, x, usecache = usecache, tax_SU = tax_SU, get_tree = get_tree), error=function(x) NULL)
  }, .progress = "text")

  #The sample_data has been corrupted by doing the merge (names get messed up and duplicated), so just regrab it with another lapply/rbind
  samp_dat <- lapply(accessions, function(x) mgnify_get_single_analysis_metadata(client, x, usecache = usecache ))
  if (returnLists){
    list(phyloseq_objects=ps_list, sample_metadata = samp_dat)
  }else{}

  #   #first of all, check to see that if we wanted them, we got trees for ALL the phyloseq objects.
  #   #If trees are present in any of the phyloseq objects during merging, then any OTUs not in a tree
  #   #(e.g. if any phyloseq objects do NOT contain a tree) will not be included in the merged output.
  #
  #   if(get_tree){
  #     if (any(is.na(lapply(ps_list, function(x) x@phy_tree)))){
  #       warning("Phylogenetic tree retrieval was requested but some of the analyses do not include phylogenetic trees. Results should be used with caution.")
  #     }
  #   }
  #
  #   #This is too slow for large datasets
  #   #full_ps <- do.call(phyloseq::merge_phyloseq, ps_list)
  #   #so:
  #   #a divide-and-conquer approach to merge_phyloseq seems to work best, hence the following
  #   #code which splits the full list into sublists and merges them seperately, then repeats until all are joined.
  #   curlist=ps_list
  #   while(length(curlist) > 1){
  #     #Lists of length 10 seem to work well
  #     sublist=split(curlist, seq_along(curlist) %/% 10)
  #     curlist <- lapply(sublist, function(x){
  #     do.call(phyloseq::merge_phyloseq,x)
  #     })
  #   }
  #
  #   #By this point curlist isn't a list, it's a phyloseq object...
  #   full_ps <- curlist[[1]]
  #   sample_metadata_df <- do.call(dplyr::bind_rows, samp_dat)
  #   rownames(sample_metadata_df) <- sample_metadata_df$analysis_accession
  #   phyloseq::sample_data(full_ps) <- sample_metadata_df
  #
  #   # col_data <- sample_data(full_ps)
  #   assay_data <- otu_table(full_ps)
  #   row_data <- tax_table(full_ps)
  #   col_data <- sample_data(full_ps)
  #   row_tree <- NULL
  #   col_tree <- NULL
  #   tse <-TreeSummarizedExperiment(assays=list(abundance=assay_data))
  #   if (!is.null(row_tree)){
  #   	rowTree(tse) <- row_tree
  #   }
  #   if (!is.null(col_tree)){
  #   	colTree(tse) <- col_tree
  #   }
  #   if (!is.null(row_data)){
  #   	rowData(tse) <- row_data
  #   }
  #   if (!is.null(col_data)){
  #     colData(tse) <- DataFrame(as.matrix(col_data))
  #   }
  #   tse
  # }
}
