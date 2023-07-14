context("getResult")
test_that("getResult", {
    # Test that input check caches wrong arguments.
    # List of arguments with correct values
    mg <- MgnifyClient()
    var <- list(
        x = list(mg),
        accession = list("test", "studies", c("studies", "assembly"),
                         "TreeSE", c("TreeSE", "phyloseq"), "taxonomy-ssu",
                         c("taxonomy-ssu", "go-slim")),
        output = list("TreeSE"),
        get.taxa = list(TRUE, FALSE),
        get.func = list(TRUE, FALSE, "test", "studies",
                        c("studies", "assembly"), "TreeSE",
                        c("TreeSE", "phyloseq"), "taxonomy-ssu",
                        c("taxonomy-ssu", "go-slim")),
        use.cache = list(TRUE, FALSE),
        verbose = list(TRUE, FALSE)
    )
    var <- .wrong_arguments(var)
    # Loop through rows, all variable pairs should end up to error
    for( i in seq_len(nrow(var)) ){
        expect_error(
            suppressWarnings(
            getResult(
                x = var[i, 1][[1]],
                accession = var[i, 2][[1]],
                output = var[i, 3][[1]],
                get.taxa = var[i, 4][[1]],
                get.func = var[i, 5][[1]],
                use.cache = var[i, 6][[1]],
                verbose = var[i, 7][[1]],
            ))
        )
    }
    # Require internet access
    skip_if(httr::http_error("https://www.ebi.ac.uk/metagenomics/api/v1"))

    # Test that only functional data is fetched based on certain accession ID.
    # Get data as list of data.frames
    res <- getResult(
        mg, "MGYA00097621",
        get.taxa = FALSE,
        output = "list",
        get.func = TRUE,
        verbose = FALSE)
    expect_true(is.list(res))
    expect_true("go-terms" %in% names(res))
    expect_true(
        is.character(res$`interpro-identifiers`$analysis) &&
        is.character(res$`interpro-identifiers`$description) &&
        is.numeric(res$`interpro-identifiers`$count))

    # Test that microbial profiling data and functional data is fetched. Get
    # data as MAE. Fetch also trees. Check that all data is is in correct place
    # and is correct.
    skip <- TRUE
    if (!skip) {
    res <- getResult(
        mg, "MGYA00097621",
        get.tree=TRUE,
        get.func=TRUE,
        verbose=FALSE)
    expect_true(class(res) == "MultiAssayExperiment")
    expect_true(class(res[[1]]) == "TreeSummarizedExperiment")
    expect_true(!is.null(rowTree(res[["microbiota"]])))
    expect_true(is.matrix(assay(res[[1]])))
    expect_true("microbiota" %in% names(res) && "go-terms" %in% names(res))
    expect_true(is.matrix(assay(res[[2]])))
    expect_true(is.matrix(assay(res[[3]])))
    expect_equal(assay(res[["go-slim"]])["GO:1990204", 1][[1]], 929)
    expect_equal(colnames(res[[1]]), colnames(res[[2]]))
    expect_equal(colnames(res[[3]]), colnames(res[[2]]))
    }
})
