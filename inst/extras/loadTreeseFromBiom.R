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
