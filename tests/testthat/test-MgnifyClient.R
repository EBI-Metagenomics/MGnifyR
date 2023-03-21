context("MgnifyClient")
test_that("MgnifyClient", {
    # Test that input check caches wrong arguments.
    # List of arguments with correct values
    mg <- MgnifyClient()
    var <- list(
        username = list("test", "study", "TreeSE", "taxonomy-ssu", NULL),
        password = list("test", "study", "TreeSE", "taxonomy-ssu", NULL),
        useCache = list(TRUE, FALSE),
        cacheDir = list("test", "study", "TreeSE", "taxonomy-ssu", NULL),
        warnings = list(TRUE, FALSE),
        useMemCache = list(TRUE, FALSE),
        url = list("test", "study", "TreeSE", "taxonomy-ssu")
    )
    var <- .wrong_arguments(var)
    # Loop through rows, all variable pairs should end up to error
    for(i in seq_len(nrow(var)) ){
        expect_error(
            MgnifyClient(
                username = var[i, 2][[1]],
                password = var[i, 3][[1]],
                useCache = var[i, 4][[1]],
                cacheDir = var[i, 5][[1]],
                warnings = var[i, 6][[1]],
                useMemCache = var[i, 7][[1]],
                url = var[i, 8][[1]],
            )
        )
    }
    # Test that slots are updated. Change arguments --> check that values
    # of slots correspond argument.
    mg <- MgnifyClient(
        useCache = TRUE,
        cacheDir = "test",
        warnings = FALSE,
        useMemCache = TRUE,
        url = "test"
    )
    expect_equal(mg@cacheDir, "test")
    expect_equal(mg@warnings, FALSE)
    expect_equal(mg@useMemCache, TRUE)
    expect_equal(mg@url, "test")
    mg <- MgnifyClient(
        useCache = FALSE,
        cacheDir = "test",
        warnings = TRUE,
        useMemCache = FALSE,
    )
    expect_true(is.na(mg@cacheDir))
    expect_equal(mg@warnings, TRUE)
    expect_equal(mg@useMemCache, FALSE)
    # Require internet access
    skip_if(httr::http_error("https://www.ebi.ac.uk/metagenomics/api/v1"))
    # Test that error occurs when wrong username/password is used in
    # authentication
    expect_error(MgnifyClient(
        username = "not_work",
        password = "not_work"
    ))
})
