context("getFile")
test_that("getFile", {
    # Test that input check caches wrong arguments.
    # List of arguments with correct values
    mg <- MgnifyClient()
    var <- list(
        x = list(mg),
        url = list("test", "studies", "TreeSE", "taxonomy-ssu"),
        file = list("test", "studies", "TreeSE", "taxonomy-ssu", NULL),
        read.func = list(NULL),
        use.cache = list(TRUE, FALSE)
    )
    var <- .wrong_arguments(var)
    # Loop through rows, all variable pairs should end up to error
    for( i in seq_len(nrow(var)) ){
        expect_error(
            getFile(
                x = var[i, 1][[1]],
                url = var[i, 2][[1]],
                file = var[i, 3][[1]],
                as.df = var[i, 4][[1]],
                read.func = var[i, 5][[1]],
                use.cache = var[i, 6][[1]]
            )
        )
    }
    # List of arguments with correct values
    mg <- MgnifyClient()
    var <- list(
        x = list(mg),
        accession = list("test", "studies", c("studies", "assembly"), "TreeSE",
                         c("TreeSE", "phyloseq"), "taxonomy-ssu",
                         c("taxonomy-ssu", "go-slim"), NULL),
        type = list("studies"),
        use.cache = list(TRUE, FALSE),
        verbose = list(TRUE, FALSE)
    )
    var <- .wrong_arguments(var)
    # Loop through rows, all variable pairs should end up to error
    for( i in seq_len(nrow(var)) ){
        expect_error(
            searchFile(
                x = var[i, 1][[1]],
                accession = var[i, 2][[1]],
                type = var[i, 3][[1]],
                use.cache = var[i, 4][[1]],
                verbose = var[i, 5][[1]]
            )
        )

    }
    # Require internet access
    skip_if(httr::http_error("https://www.ebi.ac.uk/metagenomics/api/v1"))

    # Test that df is returned even if accession ID is not correct
    res <- searchFile(mg, type="assemblies", accession="random", verbose = FALSE)
    expect_true(is.data.frame(res))

    # Test that file search is done correctly based on accession ID.
    # Use studies as type
    mg <- MgnifyClient()
    res <- searchFile(
        mg, type="studies", accession="MGYS00005292", verbose = FALSE)
    expect_true(all(res$type == "studies"))
    expect_true(is.data.frame(res))
    expect_true(grepl("https", res$download_url[1]))

    # Test that correct file is fetched based on provided url.
    res <- getFile(mg, res$download_url[1])
    # Result is stored in a path which is returned
    expect_true(file.exists(res))
})
