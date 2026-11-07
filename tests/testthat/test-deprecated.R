context("deprecated")
# The functions of earlier versions are deprecated, but they still work and
# tell what to use instead. The API is mocked, so these tests do not use
# MGnify.

test_that("MgnifyClient is not needed anymore", {
    old <- options(MGnifyR.cache = NULL)
    on.exit(options(old))
    expect_warning(
        res <- MgnifyClient(useCache = TRUE, cacheDir = "my_cache"),
        "mgnifyLogin")
    expect_null(res)
    # The cache directory is set with the option
    expect_equal(getOption("MGnifyR.cache"), "my_cache")
    # Accessors of the client are defunct
    expect_error(databaseUrl(), "defunct")
    expect_error(cacheDir(), "options\\(MGnifyR.cache\\)")
    expect_error(`useCache<-`(NULL, TRUE), "defunct")
    expect_error(mgnify_client(), "MgnifyClient")
    expect_error(mgnify_get_analyses_results(), "fetchAnalyses")
})

test_that("getResult uses fetchAnalyses", {
    type_used <- character(0)
    local_mocked_bindings(fetchAnalyses = function(analysis.ids, type, ...){
        type_used <<- c(type_used, type)
        analysis.ids
    })
    expect_warning(res <- getResult(NULL, "MGYA1", get.func = FALSE),
        "deprecated")
    expect_equal(res, "MGYA1")
    # Old taxonomy names, e.g., "taxonomy-lsu", are the types
    suppressWarnings(getResult(NULL, "MGYA1", get.taxa = "taxonomy-lsu"))
    # Unknown taxonomy gives SSU
    suppressWarnings(getResult(NULL, "MGYA1", get.taxa = "taxonomy"))
    expect_equal(type_used, c("SSU", "lsu", "SSU"))
    # Functional data is not returned as TreeSE anymore
    expect_warning(
        getResult(NULL, "MGYA1", get.func = "go-slim"),
        "Only taxonomic profiles")
})

test_that("searchAnalysis and getMetadata use searchAnalyses", {
    analyses <- data.frame(
        analysis_accession = c("MGYA1", "MGYA2", "MGYA3"),
        study_accession = c("MGYS1", "MGYS1", "MGYS2"),
        sample_accession = c("SAMEA1", "SAMEA2", "SAMEA1"))
    local_mocked_bindings(
        searchAnalyses = function(study.ids, ...){
            analyses[analyses[["study_accession"]] %in% study.ids, ]
        },
        # Sample given by its ENA accession
        .v2_get_json = function(urls) list(list(
            accession = "SAMEA1", studies = list(list(accession = "MGYS1"),
                list(accession = "MGYS2")))),
        .v2_analyses_by_id = function(ids, ...){
            # The old function returned all the metadata
            expect_true(all(unlist(list(...))))
            analyses[analyses[["analysis_accession"]] %in% ids, ]
        })
    # Analyses of studies are named by the study
    expect_warning(res <- searchAnalysis(NULL, "studies", "MGYS1"),
        "searchAnalyses")
    expect_equal(res, c(MGYS1 = "MGYA1", MGYS1 = "MGYA2"))
    # Analyses of samples are named by the given sample accession
    expect_warning(res <- searchAnalysis(NULL, "samples", "ERS1"))
    expect_equal(res, c(ERS1 = "MGYA1", ERS1 = "MGYA3"))
    expect_error(
        suppressWarnings(searchAnalysis(NULL, "runs", "ERR1")),
        "'type' must be")
    expect_warning(res <- getMetadata(NULL, c("MGYA2", "MGYA2")),
        "searchAnalyses")
    expect_equal(rownames(res), "MGYA2")
})

test_that("doQuery and getData query any endpoint", {
    local_mocked_bindings(
        # Detail endpoints have a single item, and list endpoints have pages
        .v2_get_json = function(urls) lapply(urls, function(url){
            if( grepl("studies/MGYS1$", url) ){
                return(list(accession = "MGYS1", title = "a"))
            }
            list(count = 2L, items = list())
        }),
        .v2_get_all = function(paths, ...) list(list(
            list(accession = "MGYS2"), list(accession = "MGYS3"))))
    expect_warning(res <- doQuery(NULL, "studies", "MGYS1"), "searchStudies")
    expect_equal(res[["title"]], "a")
    expect_warning(res <- doQuery(NULL, "studies", as.df = FALSE))
    expect_length(res, 2L)
    expect_warning(
        res <- getData(NULL, "analyses", "studies", "MGYS1"),
        "searchAnalyses")
    expect_equal(res[["accession"]], c("MGYS2", "MGYS3"))
})

test_that("searchFile and getFile list and download files", {
    local_mocked_bindings(.v2_get_json = function(urls) list(
        list(accession = "MGYS1", downloads = list(
            list(alias = "a.tsv", url = "https://x/a.tsv"))),
        NULL))
    expect_message(expect_warning(
        res <- searchFile(NULL, c("MGYS1", "MGYS2"), type = "studies"),
        "deprecated"), "MGYS2")
    expect_equal(res[["accession"]], "MGYS1")
    expect_equal(res[["url"]], "https://x/a.tsv")

    path <- tempfile(fileext = ".tsv")
    writeLines("x", path)
    local_mocked_bindings(.v2_download = function(urls, ext){
        expect_equal(ext, ".tsv")
        ifelse(grepl("missing", urls), NA_character_, path)
    })
    expect_warning(res <- getFile(NULL, "https://x/a.tsv"), "deprecated")
    expect_equal(res, path)
    copy <- tempfile()
    suppressWarnings(getFile(NULL, "https://x/a.tsv", file = copy))
    expect_equal(readLines(copy), "x")
    res <- suppressWarnings(
        getFile(NULL, "https://x/a.tsv", read.func = readLines))
    expect_equal(res, "x")
    expect_error(
        suppressWarnings(getFile(NULL, "https://x/missing.tsv")),
        "Failed to download")
})
