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
        url = "test"
    )
    expect_equal(cacheDir(mg), "test/.MGnifyR_cache")
    expect_equal(showWarnings(mg), FALSE)
    expect_equal(databaseUrl(mg), "test")
    mg <- MgnifyClient(
        useCache = FALSE,
        cacheDir = "test",
        showWarnings = TRUE
    )
    expect_true(!is.na(cacheDir(mg)))
    expect_equal(showWarnings(mg), TRUE)
    # Require internet access
    skip_if(httr::http_error("https://www.ebi.ac.uk/metagenomics/api/v1"))
    # Test that error occurs when wrong username/password is used in
    # authentication
    expect_error(MgnifyClient(username = "not_work", password = "not_work"))
    expect_error(
        MgnifyClient(
            username = "not_work", password = "not_work", url = "not_work"))
})
