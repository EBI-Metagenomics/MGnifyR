loadTreeseFromBiom <- function (BIOMfilename, treefilename = NULL, refseqfilename = NULL,
          refseqFunction = readDNAStringSet, refseqArgs = NULL, parseFunction = parse_taxonomy_default,
          parallel = FALSE, version = 1, ...)
{
  argumentlist <- list()
  if (class(BIOMfilename) == "character") {
    x = biomformat::read_biom(biom_file = BIOMfilename)
  }
  else{
    if (class(BIOMfilename) == "biom") {
     x = BIOMfilename
    } else {
       stop("import_biom requires a 'character' string to a biom file or a 'biom-class' object")
    }
  }

  otutab = otu_table(as(biomformat::biom_data(x), "matrix"), taxa_are_rows = TRUE)
  argumentlist <- c(argumentlist, list(otutab))
  if (all(sapply(sapply(x$rows, function(i) {
    i$metadata
  }), is.null))) {
    taxtab <- NULL
  }
  else {
    taxlist = lapply(x$rows, function(i) {
      parseFunction(i$metadata$taxonomy)
    })
    names(taxlist) = sapply(x$rows, function(i) {
      i$id
    })
    taxtab = build_tax_table(taxlist)
  }
  argumentlist <- c(argumentlist, list(taxtab))
  if (is.null(biomformat::sample_metadata(x))) {
    samdata <- NULL
  }
  else {
    samdata = sample_data(biomformat::sample_metadata(x))
  }
  argumentlist <- c(argumentlist, list(samdata))
  tree <- NULL
  if (!is.null(treefilename)) {
    if (inherits(treefilename, "phylo")) {
      tree = treefilename
    }
    else {
      tree <- read_tree(treefilename, ...)
    }
    if (is.null(tree)) {
      warning("treefilename failed import. It not included.")
    }
    else {
      argumentlist <- c(argumentlist, list(tree))
    }
  }

  assay_data <- otutab
  row_data <- taxtab
  col_data <- samdata
  row_tree <- tree
  tse <-TreeSummarizedExperiment(assays=list(abundance=assay_data))
  if (!is.null(row_tree)){
    rowTree(tse) <- row_tree
  }
  if (!is.null(row_data)){
    rowData(tse) <- row_data
  }
  if (!is.null(col_data)){
    colData(tse) <- DataFrame(as.matrix(col_data))
  }
  tse
}

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

  #Load in a TreeSummarizedExperiment object
  tse <- loadTreeseFromBiom(biom_path)

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
