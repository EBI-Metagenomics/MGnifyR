context("getData")
test_that("getData", {
    # Test that input check caches wrong arguments.
    mg <- MgnifyClient(useCache = FALSE)

    expect_error(getData(1))
    expect_error(getData("test"))
    expect_error(getData(TRUE))

    expect_error(getData(mg, type = 1))
    expect_error(getData(mg, type = TRUE))
    expect_error(getData(mg, type = NULL))
    expect_error(getData(mg, type = c("type", "type")))

    expect_error(getData(
        mg, type = "kegg-modules", accession.type = "analyses", accesion = 1))
    expect_error(getData(
        mg, type = "kegg-modules", accession.type = "analyses",
        accesion = TRUE))
    expect_error(getData(
        mg, type = "kegg-modules", accession.type = "analyses",
        accesion = NULL))

    expect_error(getData(
        mg, type = "kegg-modules", accession = c("MGYA00642773"),
        accesion.type = 1))
    expect_error(getData(
        mg, type = "kegg-modules",
        accession = c("MGYA00642773", "MGYA00642774"), accesion.type = TRUE))
    expect_error(getData(
        mg, type = "kegg-modules",
        accession = c("MGYA00642773", "MGYA00642774"), accesion.type = NULL))
    expect_error(getData(
        mg, type = "kegg-modules",
        accession = c("MGYA00642773", "MGYA00642774"), accesion.type = c("type", "type")))

    expect_error(getData(
        mg, type = "kegg-modules", accession = c("MGYA00642773"),
        accesion.type = c("type"), as.df = "test"))
    expect_error(getData(
        mg, type = "kegg-modules", accession = c("MGYA00642773"),
        accesion.type = c("type"), as.df = 1))
    expect_error(getData(
        mg, type = "kegg-modules", accession = c("MGYA00642773"),
        accesion.type = c("type"), as.df = c(TRUE, TRUE)))
    expect_error(getData(
        mg, type = "kegg-modules", accession = c("MGYA00642773"),
        accesion.type = c("type"), as.df = NULL))

    # Require internet access
    skip_if(httr::http_error("https://www.ebi.ac.uk/metagenomics/api/v1"))

    type <- "kegg-modules"
    res <- getData(
      mg, type = type, accession = "MGYA00642773", accession.type = "analyses")
    expect_is(res, "data.frame")
    expect_true( all(res[["type"]] == type) )
})
