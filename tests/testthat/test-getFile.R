context("getFile")
test_that("getFile", {
    # Test that input check caches wrong arguments.
    mg <- MgnifyClient()

    expect_error(getFile(10))
    expect_error(getFile(TRUE))
    expect_error(getFile(NULL))

    expect_error(getFile(mg, url = 10))
    expect_error(getFile(mg, url = TRUE))
    expect_error(getFile(mg, url = c("test", "test")))

    expect_error(getFile(mg, url = "test", read.func = 10))
    expect_error(getFile(mg, url = "test", read.func = TRUE))

    expect_error(getFile(mg, url = "test", use.cache = 10))
    expect_error(getFile(mg, url = "test", use.cache = TRUE))
    expect_error(getFile(mg, url = "test", use.cache = c("test", "test")))

    expect_error(getFile(mg, url = "taxonomy--ssu", use.cache = 10))
    expect_error(getFile(mg, url = "test", use.cache = test))

    expect_error(searchFile(10))
    expect_error(searchFile(TRUE))
    expect_error(searchFile(NULL))

    expect_error(searchFile(mg, accession = TRUE))
    expect_error(searchFile(mg, accession = 1))
    expect_error(searchFile(mg, accession = NULL))

    expect_error(searchFile(mg, accession = "test", type = 1))
    expect_error(searchFile(mg, accession = "test", type = TRUE))
    expect_error(searchFile(mg, accession = "test", c("samples", "analyses")))

    expect_error(searchFile(mg, accession = "test", type = "samples", use.cache = NULL))
    expect_error(searchFile(mg, accession = "test", type = "samples", use.cache = 1))
    expect_error(searchFile( mg, accession = "test", type = "samples", use.cache = c(TRUE, FALSE)))

    expect_error(searchFile(mg, accession = "test", type = "samples", verbose = NULL))
    expect_error(searchFile(mg, accession = "test", type = "samples", verbose = 1))
    expect_error(searchFile( mg, accession = "test", type = "samples", verbose = c(TRUE, FALSE)))

    # Require internet access
    skip_if(httr::http_error("https://www.ebi.ac.uk/metagenomics/api/v1"))

    # Expect error because url is incorrect
    expect_error(getFile(mg, url = "test"))

    # Test that df is returned even if accession ID is not correct
    expect_warning(
    res <- searchFile(mg, type = "assemblies", accession = "random")
    )
    expect_true(is.data.frame(res))

    # Test that file search is done correctly based on accession ID.
    # Use studies as type
    res <- searchFile(mg, type = "studies", accession = "MGYS00005292", verbose = FALSE)
    expect_true(all(res$type == "studies"))
    expect_true(is.data.frame(res))
    expect_true(grepl("https", res$download_url[1]))

    # Test that correct file is fetched based on provided url.
    res <- getFile(mg, res$download_url[1])
    # Result is stored in a path which is returned
    expect_true(file.exists(res))
})
