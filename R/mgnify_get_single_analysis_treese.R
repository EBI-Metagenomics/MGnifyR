# get a single biom file and convert it to a TreeSummarizedExperiment object
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
    warning("Unable to locate requested taxonomy type ",tax_SU,". This is likely due to the current analysis having been performed on an older version of the MGnify pipeline. The available BIOM file will be used instead.")
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

  #Load in the TreeSummarizedExperiment object
  tse <- loadTreeseFromBiom(biom_path)
  # tse <- loadFromBiom(biom_path)

  #Need to check if the taxonomy was parsed correctly - depending on the pipeline it may need a bit of help:
  mia::checkTaxonomy(tse)

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
