context("searchAnalysis")
test_that("searchAnalysis", {
    # Test that input check caches wrong arguments.
    mg <- MgnifyClient()

    expect_error(searchAnalysis(TRUE))
    expect_error(searchAnalysis("test"))
    expect_error(searchAnalysis(NULL))
    expect_error(searchAnalysis(1))

    expect_error(searchAnalysis(mg, type = TRUE, accession = "test"))
    expect_error(searchAnalysis(mg, type = "test", accession = "test"))
    expect_error(searchAnalysis(mg, type = NULL, accession = "test"))
    expect_error(searchAnalysis(mg, type = c("studies", "samples", accession = "test")))

    expect_error(searchAnalysis(mg, type = "studies", accession = TRUE))
    expect_error(searchAnalysis(mg, type = "studies", accession = NULL))
    expect_error(searchAnalysis(mg, type = "studies", accession = 1))
    expect_error(searchAnalysis(mg, type = "studies", accession = c(TRUE, FALSE)))

    expect_error(searchAnalysis(mg, type = "studies", accession = "test", use.cache = NULL))
    expect_error(searchAnalysis(mg, type = "studies", accession = "test", use.cache = 1))
    expect_error(searchAnalysis(mg, type = "studies", accession = "test", use.cache = "test"))

    expect_error(searchAnalysis(mg, type = "studies", accession = "test", verbose = NULL))
    expect_error(searchAnalysis(mg, type = "studies", accession = "test", verbose = 1))
    expect_error(searchAnalysis(mg, type = "studies", accession = "test", verbose = "test"))

    # Require internet access
    skip_if(httr::http_error("https://www.ebi.ac.uk/metagenomics/api/v1"))

    # Test that correct analysis IDs are found based on study ID.
    expect_warning(res <- searchAnalysis(mg, "studies", "MGYA00097621", verbose = FALSE))
    expect_true(is.null(res))
    res <- searchAnalysis(mg, "studies", "MGYS00005058", verbose = FALSE)
    expect_true(length(res) > 0)
    expect_true("MGYA00377528" %in% res)
    # Test that correct analysis IDs are found based on sample ID.
    expect_warning(res <- searchAnalysis(mg, "samples", "MGYA00097621", verbose = FALSE))
    expect_true(is.null(res))
    res <-  searchAnalysis(mg, "samples", "ERS2161777", verbose = FALSE)
    expect_true(length(res) > 0)
    expect_true("MGYA00293854" %in% res)
})
