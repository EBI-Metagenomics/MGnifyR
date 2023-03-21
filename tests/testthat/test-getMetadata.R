context("getMetadata")
test_that("getMetadata", {
    # Test that input check caches wrong arguments.
    # List of arguments with correct values
    mg <- MgnifyClient()
    var <- list(
        x = list(mg),
        accession = list("test", "studies", c("studies", "assembly"),
                         "TreeSE", c("TreeSE", "phyloseq"), "taxonomy-ssu",
                         c("taxonomy-ssu", "go-slim")),
        use.cache = list(TRUE, FALSE),
        verbose = list(TRUE, FALSE)
    )
    var <- .wrong_arguments(var)
    # Loop through rows, all variable pairs should end up to error
    for( i in seq_len(nrow(var)) ){
        expect_error(
            getMetadata(
                x = var[i, 1][[1]],
                accession = var[i, 2][[1]],
                use.cache = var[i, 3][[1]],
                verbose = var[i, 4][[1]]
            )
        )
    }
    # Require internet access
    skip_if(httr::http_error("https://www.ebi.ac.uk/metagenomics/api/v1"))

    # Test that correct metadata is fetched based on certain accession ID.
    res <- getMetadata(mg, "MGYA00097621", verbose = FALSE)
    expect_equal(nrow(res), 1)
    expect_true(ncol(res) > 1)
    expect_equal(rownames(res)[1] , "MGYA00097621")
    expect_equal(res$run_accession, "ERR1811651")
    # When metadata is not found, should give a warning and the result should
    # be empty tibble
    expect_warning(
    res <- getMetadata(mg, "MGYS00005292", verbose = FALSE)
    )
    expect_true(ncol(res) == 0 && nrow(res) == 0)
})
