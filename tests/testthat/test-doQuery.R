context("doQuery")
test_that("doQuery", {
    # Test that input check caches wrong arguments.
    # List of arguments with correct values
    mg <- MgnifyClient()
    var <- list(
        x = list(mg),
        type = list("studies"),
        accession = list("test", "studies", c("studies", "assembly"), "TreeSE",
                         c("TreeSE", "phyloseq"), "taxonomy-ssu",
                         c("taxonomy-ssu", "go-slim"), NULL),
        as.df = list(TRUE, FALSE),
        max.hits = list(NULL, 0, 1, 16),
        use.cache = list(TRUE, FALSE)
    )
    var <- .wrong_arguments(var)
    # Loop through rows, all variable pairs should end up to error
    for( i in seq_len(nrow(var)) ){
        expect_error(
            doQuery(
                x = var[i, 1][[1]],
                type = var[i, 2][[1]],
                accession = var[i, 3][[1]],
                as.df = var[i, 4][[1]],
                max.hits = var[i, 5][[1]],
                use.cache = var[i, 6][[1]]
            )
        )
    }
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
    expect_equal(query2$bioproject, query$MGYS00005292$attributes$bioproject)
})
