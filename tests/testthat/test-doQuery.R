context("doQuery")
test_that("doQuery", {
    # Test that input check caches wrong arguments.
    mg <- MgnifyClient()

    # Expect errors when input is wrong
    expect_error(doQuery("test"))
    expect_error(doQuery(TRUE))
    expect_error(doQuery(1))

    expect_error(doQuery(mg, type = 1))
    expect_error(doQuery(mg, type = "test"))
    expect_error(doQuery(mg, type = TRUE))
    expect_error(doQuery(mg, type = c("studies", "samples")))

    expect_error(doQuery(mg, type = "studies", accession = 1))
    expect_error(doQuery(mg, type = "studies", accession = TRUE))
    expect_error(doQuery(mg, type = "studies", accession = c(1, 2)))

    expect_error(doQuery(mg, type = "studies", accession = "test", as.df = NULL))
    expect_error(doQuery(mg, type = "studies", accession = "test", as.df = 1))
    expect_error(doQuery(mg, type = "studies", accession = "test", as.df = c(TRUE, FALSE)))

    expect_error(doQuery(mg, type = "studies", accession = "test", max.hits = TRUE))
    expect_error(doQuery(mg, type = "studies", accession = "test", max.hits = -100))
    expect_error(doQuery(mg, type = "studies", accession = "test", max.hits = c(1, 2)))
    expect_error(doQuery(mg, type = "studies", accession = "test", max.hits = 1.5))

    expect_error(doQuery(mg, type = "studies", accession = "test", use.cache = NULL))
    expect_error(doQuery(mg, type = "studies", accession = "test", use.cache = 1))
    expect_error(doQuery(mg, type = "studies", accession = "test", use.cache = c(TRUE, FALSE)))

    # Require internet access
    skip_if(httr::http_error("https://www.ebi.ac.uk/metagenomics/api/v1"))
    # Test that studies are searched based on certain accession ID, get result
    # as list, choose max hits
    query <- doQuery(mg, "studies", "MGYS00005292", max.hits = 1, as.df = FALSE)
    expect_true(is.list(query))
    expect_true(names(query) %in% "MGYS00005292")
    expect_true(query$MGYS00005292$type == "studies")
    # Test that runs are searched, get result as df, choose max hits
    query2 <- doQuery(mg, "studies", "MGYS00005292", max.hits = 1)
    expect_true(is.data.frame(query2))
    expect_equal(query2$bioproject,
                 query$MGYS00005292$attributes$bioproject)
})
