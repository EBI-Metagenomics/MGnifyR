#' @importFrom urltools parameters
#' @importFrom httr GET
#' @importFrom httr write_disk
#' @importFrom phyloseq import_biom
#' @importFrom phyloseq tax_table
#' @importFrom phyloseq parse_taxonomy_greengenes
#' @importFrom phyloseq sample_names
#' @importFrom phyloseq sample_data
#' @importFrom phyloseq phy_tree
#' @importFrom ape read.tree
#' @importFrom utils tail

#UPDATE ME SO THAT TREES (if available) GET GRABBED AS WELL!!!
# Not exported - get a single biom file and convert it to a phyloseq object.
mgnify_get_single_analysis_phyloseq <- function(client=NULL, accession, usecache=T, downloadDIR=NULL, tax_SU="SSU", get_tree=FALSE){
    metadata_df <- mgnify_get_single_analysis_metadata(client, accession, usecache=usecache)

    analysis_data <- mgnify_retrieve_json(client, paste("analyses",accession,sep="/"), usecache = usecache)
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
        downloadDIR <- paste(client@cache_dir,"biom_files",sep="/")
        dir.create(downloadDIR, recursive = T, showWarnings = client@warnings)
    }
    #Clear out any ?params after the main location - don't need them for this
    urltools::parameters(biom_url) <- NULL

    fname <- tail(strsplit(biom_url, '/')[[1]], n=1)
    biom_path <- paste(downloadDIR, fname, sep="/")

    ## Quick check to see if we should clear the disk cache ~for this specific call~ - used for debugging
    # and when MGnify breaks
    if(usecache & client@clear_cache){
        message(paste("clear_cache is TRUE: deleting ",biom_path, sep=""))
        tryCatch(unlink(biom_path), error=warning)
    }

    if (! file.exists(biom_path)){#} | !use_downloads ){
        httr::GET(biom_url, httr::write_disk(biom_path, overwrite = T))
    }
    #Load in the phlyloseq object
    psobj <- phyloseq::import_biom(biom_path)
    #Need to check if the taxonomy was parsed correctly - depending on the pipeline it may need a bit of help:
    if (ncol(phyloseq::tax_table(psobj)) == 1){
        psobj <- phyloseq::import_biom(biom_path, parseFunction = phyloseq::parse_taxonomy_qiime)
    }
    if(! "Kingdom" %in% colnames(phyloseq::tax_table(psobj))){
        psobj <- phyloseq::import_biom(biom_path, parseFunction = phyloseq::parse_taxonomy_greengenes)
    }

    #The biom files have a single column of unknown name - maybe it's the original sample name?.
    # It's rewritten as sample_run_analysis accession, with
    # the original value stored in the sample_data (just in case it's needed later)
    orig_samp_name <- phyloseq::sample_names(psobj)[[1]]
    newsampname <- rownames(metadata_df)[1]
    metadata_df[1,"orig_samp_name"] <- orig_samp_name
    phyloseq::sample_names(psobj) <- newsampname
    phyloseq::sample_data(psobj) <- metadata_df

    #Finally, do we want to the phylogenetic tree? If so, is it there?
    if(get_tree){
        #is there a tree?
        tvec <- grepl('Phylogenetic tree', sapply(analysis_downloads, function(x) x$attributes$`description`$label))
        if(any(tvec)){
            tree_url <- analysis_downloads[tvec][[1]]$links$self
            #Clear out any ?params after the main location - don't need them for this
            urltools::parameters(tree_url) <- NULL
            #@importFrom utils tail
            fname=tail(strsplit(tree_url, '/')[[1]], n=1)
            tree_path = paste(downloadDIR, fname, sep="/")

            ## Quick check to see if we should clear the disk cache ~for this specific call~ - used for debugging
            # and when MGnify breaks
            if(usecache & client@clear_cache){
                message(paste("clear_cache is TRUE: deleting ",tree_path, sep=""))
                tryCatch(unlink(tree_path), error=warning)
            }

            if (! file.exists(tree_path)){#} | !use_downloads ){
                httr::GET(tree_url, httr::write_disk(tree_path, overwrite = T ))
            }
        }
        phylo_tree <- ape::read.tree(tree_path)
        phyloseq::phy_tree(psobj) <- phylo_tree
    }
    psobj
}
