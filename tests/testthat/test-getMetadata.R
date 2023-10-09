context("getMetadata")
test_that("getMetadata", {
    # Test that input check caches wrong arguments.
    mg <- MgnifyClient()

    expect_error(getMetadata(1))
    expect_error(getMetadata("test"))
    expect_error(getMetadata(TRUE))

    expect_error(getMetadata(mg, accession = TRUE))
    expect_error(getMetadata(mg, accession = 1))
    expect_error(getMetadata(mg, accession = c(1, 2)))

    expect_error(getMetadata(mg, accession = "test", use.cache = NULL))
    expect_error(getMetadata(mg, accession = "test", use.cache = 1))
    expect_error(getMetadata(mg, accession = "test", use.cache = c(TRUE, FALSE)))

    expect_error(getMetadata(mg, accession = "test", verbose = NULL))
    expect_error(getMetadata(mg, accession = "test", verbose = 1))
    expect_error(getMetadata(mg, accession = "test", verbose = c(TRUE, FALSE)))

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
    expect_warning(res <- getMetadata(mg, "MGYS00005292", verbose = FALSE))
    expect_true(ncol(res) == 0 && nrow(res) == 0)
})
