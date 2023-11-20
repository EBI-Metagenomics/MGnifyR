context("MgnifyClient")
test_that("MgnifyClient", {
    # Test that input check caches wrong arguments.
    mg <- MgnifyClient()

    # Expect errors when input is wrong
    expect_error(MgnifyClient(useCache = 1))
    expect_error(MgnifyClient(useCache = "TRUE"))
    expect_error(MgnifyClient(useCache = c(TRUE, TRUE)))
    
    expect_error(MgnifyClient(verbose = 1))
    expect_error(MgnifyClient(verbose = "TRUE"))
    expect_error(MgnifyClient(verbose = c(TRUE, TRUE)))
    
    expect_error(MgnifyClient(showWarnings = 1))
    expect_error(MgnifyClient(showWarnings = "TRUE"))
    expect_error(MgnifyClient(showWarnings = c(TRUE, TRUE)))
    
    expect_error(MgnifyClient(clearCache = 1))
    expect_error(MgnifyClient(clearCache = "TRUE"))
    expect_error(MgnifyClient(clearCache = c(TRUE, TRUE)))

    expect_error(MgnifyClient(showWarnings = 1))
    expect_error(MgnifyClient(showWarnings = "TRUE"))
    expect_error(MgnifyClient(showWarnings = c(TRUE, TRUE)))

    expect_error(MgnifyClient(useMemCache = 1))
    expect_error(MgnifyClient(useMemCache = "TRUE"))
    expect_error(MgnifyClient(useMemCache = c(TRUE, TRUE)))

    expect_error(MgnifyClient(url = 1))
    expect_error(MgnifyClient(url = TRUE))
    expect_error(MgnifyClient(url = c("url", "url")))

    expect_error(MgnifyClient(username = 1))
    expect_error(MgnifyClient(username = TRUE))
    expect_error(MgnifyClient(username = c("url", "url")))

    expect_error(MgnifyClient(password = 1))
    expect_error(MgnifyClient(password = TRUE))
    expect_error(MgnifyClient(password = c("url", "url")))

    # Test that slots are updated. Change arguments --> check that values
    # of slots correspond argument.
    mg <- MgnifyClient(
        useCache = TRUE,
        cacheDir = "test",
        showWarnings = FALSE,
        useMemCache = TRUE,
        url = "test"
    )
    expect_equal(mg@cacheDir, "test")
    expect_equal(mg@showWarnings, FALSE)
    expect_equal(mg@useMemCache, TRUE)
    expect_equal(mg@databaseUrl, "test")
    mg <- MgnifyClient(
        useCache = FALSE,
        cacheDir = "test",
        showWarnings = TRUE,
        useMemCache = FALSE,
    )
    expect_true(!is.na(mg@cacheDir))
    expect_equal(mg@showWarnings, TRUE)
    expect_equal(mg@useMemCache, FALSE)
    # Require internet access
    skip_if(httr::http_error("https://www.ebi.ac.uk/metagenomics/api/v1"))
    # Test that error occurs when wrong username/password is used in
    # authentication
    expect_error(MgnifyClient(username = "not_work", password = "not_work"))
    expect_error(MgnifyClient(username = "not_work", password = "not_work", url = "not_work"))
})
