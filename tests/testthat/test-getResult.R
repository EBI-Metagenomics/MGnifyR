context("getResult")
test_that("getResult", {
    # Test that input check caches wrong arguments.
    mg <- MgnifyClient()

    expect_error(getResult(1))
    expect_error(getResult("test"))
    expect_error(getResult(TRUE))

    expect_error(getResult(mg, accesion = 1))
    expect_error(getResult(mg, accesion = TRUE))
    expect_error(getResult(mg, accesion = NULL))

    expect_error(getResult(mg, accession = "test", output = "test"))
    expect_error(getResult(mg, accession = "test", output = TRUE))
    expect_error(getResult(mg, accession = "test", output = 1))
    expect_error(getResult(mg, accession = "test", output = c("TreeSE", "phyloseq")))
    expect_error(getResult(mg, accession = "test", output = NULL))

    expect_error(getResult(mg, accession = "test", get.taxa = NULL))
    expect_error(getResult(mg, accession = "test", get.taxa = 1))
    expect_error(getResult(mg, accession = "test", get.taxa = c(TRUE, TRUE)))
    expect_error(getResult(mg, accession = "test", get.taxa = "test"))

    expect_error(getResult(mg, accession = "MGYA00097621", get.tree = NULL))
    expect_error(getResult(mg, accession = "MGYA00097621", get.tree = 1))
    expect_error(getResult(mg, accession = "MGYA00097621", get.tree = c(TRUE, TRUE)))
    expect_error(getResult(mg, accession = "MGYA00097621", get.tree = "test"))

    expect_error(getResult(mg, accession = "test", get.func = NULL))
    expect_error(getResult(mg, accession = "test", get.func = 1))
    expect_error(getResult(mg, accession = "test", get.func = c(TRUE, TRUE)))
    expect_error(getResult(mg, accession = "test", get.func = "test"))

    expect_error(getResult(mg, accession = "test", verbose = NULL))
    expect_error(getResult(mg, accession = "test", verbose = 1))
    expect_error(getResult(mg, accession = "test", verbose = c(TRUE, TRUE)))
    expect_error(getResult(mg, accession = "test", verbose = "test"))

    expect_error(getResult(mg, accession = "test", use.cache = NULL))
    expect_error(getResult(mg, accession = "test", use.cache = 1))
    expect_error(getResult(mg, accession = "test", use.cache = c(TRUE, TRUE)))
    expect_error(getResult(mg, accession = "test", use.cache = "test"))

    # Require internet access
    skip_if(httr::http_error("https://www.ebi.ac.uk/metagenomics/api/v1"))

    # Test that only functional data is fetched based on certain accession ID.
    # Get data as list of data.frames
    res <- getResult( mg, "MGYA00097621", get.taxa = FALSE, output = "list", get.func = TRUE, verbose = FALSE)
    expect_true(is.list(res))
    expect_true("go-terms" %in% names(res))
    expect_true(is.character(res$`interpro-identifiers`$analysis) &&
                    is.character(res$`interpro-identifiers`$description) &&
                    is.numeric(res$`interpro-identifiers`$count))

    # Test that microbial profiling data and functional data is fetched. Get
    # data as MAE. Fetch also trees. Check that all data is is in correct place
    # and is correct.
    res <- getResult(mg, "MGYA00097621", get.tree = TRUE, get.func = TRUE, verbose = FALSE)
    expect_true(class(res) == "MultiAssayExperiment")
    expect_true(class(res[[1]]) == "TreeSummarizedExperiment")
    expect_true(!is.null(rowTree(res[["microbiota"]])))
    expect_true(is.matrix(assay(res[[1]])))
    expect_true("microbiota" %in% names(res) &&
                    "go-terms" %in% names(res))
    expect_true(is.matrix(assay(res[[2]])))
    expect_true(is.matrix(assay(res[[3]])))
    expect_equal(assay(res[["go-slim"]])["GO:1990204", 1][[1]], 929)
    expect_equal(colnames(res[[1]]), colnames(res[[2]]))
    expect_equal(colnames(res[[3]]), colnames(res[[2]]))
})
