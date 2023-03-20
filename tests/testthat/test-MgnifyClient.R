context("MgnifyClient")
test_that("MgnifyClient", {
    # Test that input check caches wrong arguments.
    # List of arguments with correct values
    mg <- MgnifyClient()
    var <- list(
        username = list("test", "study", "TreeSE", "taxonomy-ssu", NULL),
        password = list("test", "study", "TreeSE", "taxonomy-ssu", NULL),
        use.cache = list(TRUE, FALSE),
        cache.dir = list("test", "study", "TreeSE", "taxonomy-ssu", NULL),
        warnings = list(TRUE, FALSE),
        use.memcache = list(TRUE, FALSE),
        url = list("test", "study", "TreeSE", "taxonomy-ssu")
    )
    var <- .wrong_arguments(var)
    # Loop through rows, all variable pairs should end up to error
    for(i in seq_len(nrow(var)) ){
        expect_error(
            MgnifyClient(
                username = var[i, 2][[1]],
                password = var[i, 3][[1]],
                use.cache = var[i, 4][[1]],
                cache.dir = var[i, 5][[1]],
                warnings = var[i, 6][[1]],
                use.memcache = var[i, 7][[1]],
                url = var[i, 8][[1]],
            )
        )
    }
    # Test that slots are updated. Change arguments --> check that values
    # of slots correspond argument.
    mg <- MgnifyClient(
        use.cache = TRUE,
        cache.dir = "test",
        warnings = FALSE,
        use.memcache = TRUE,
        url = "test"
    )
    expect_equal(mg@cacheDir, "test")
    expect_equal(mg@warnings, FALSE)
    expect_equal(mg@useMemCache, TRUE)
    expect_equal(mg@url, "test")
    mg <- MgnifyClient(
        use.cache = FALSE,
        cache.dir = "test",
        warnings = TRUE,
        use.memcache = FALSE,
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
