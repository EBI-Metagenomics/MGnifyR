context("v2_search")
test_that("Input check", {
    expect_error(searchStudies(biome = TRUE))
    expect_error(searchStudies(biome = c("root", "root")))
    expect_error(searchStudies(biome = ""))
    expect_error(searchStudies(pipeline.version = 5))
    expect_error(searchStudies(search = NA_character_))
    expect_error(searchStudies(order = c("accession", "updated_at")))
    expect_error(searchStudies(max.hits = 0))
    expect_error(searchStudies(max.hits = 1.5))
    expect_error(searchStudies(max.hits = "10"))

    expect_error(searchSamples(study.ids = 1))
    expect_error(searchSamples(study.ids = character(0)))
    expect_error(searchSamples(study.ids = data.frame(x = "MGYS1")))
    expect_error(searchSamples("MGYS00005116", biome = "root"))
    expect_error(searchSamples("MGYS00005116", search = "a"))
    expect_error(searchSamples("MGYS00005116", order = "updated_at"))

    expect_error(searchAnalyses())
    expect_error(searchAnalyses(NULL))
    expect_error(searchAnalyses(c("MGYS1", NA)), "character vector")
    expect_error(searchAnalyses("MGYS1", qc = NA), "'qc' must be")
    expect_error(
        searchAnalyses("MGYS1", biosamples = "yes"), "'biosamples' must be")

    expect_error(mgnifyLogin("user"))
    expect_error(mgnifyLogin(1, "password"))
})

test_that("Helpers", {
    expect_equal(
        .v2_norm_version(c("V5", "v5.0", "5.0", "V4.1", "6", "V6.1")),
        c("5", "5", "5", "4.1", "6", "6.1"))
    expect_equal(
        .v2_url(c("studies/", "samples/"), list(page = 2, search = NULL)),
        paste0(.V2_URL, c("studies/?page=2", "samples/?page=2")))
    expect_equal(
        .v2_url("studies/", list(biome_lineage = "root:A b")),
        paste0(.V2_URL, "studies/?biome_lineage=root%3AA%20b"))
    expect_equal(.v2_url("studies/MGYS1"), paste0(.V2_URL, "studies/MGYS1"))
    expect_equal(.v2_download(character(0), ".json"), character(0))
    expect_equal(.v2_ids(c("A", "A", "B"), "x", "accession"), c("A", "B"))
    expect_equal(.v2_ids(data.frame(accession = c("A", "A", "B")), "x",
        "accession"), c("A", "B"))
    # Nested objects become columns named by their path, arrays are JSON,
    # and file lists are left out
    item <- list(
        accession = "MGYA1", ena_accessions = list("ERR1", "SAM1"),
        run = list(accession = "ERR1", instrument_model = "x"),
        sample = list(sample_title = "title"), assembly = NULL,
        downloads = list(list(url = "a"), list(url = "b")),
        publications = list(list(pubmed_id = 1L), list(pubmed_id = 2L)),
        empty = list(nested = list(list(a = 1))))
    res <- .v2_flatten(item)
    expect_equal(res[["accession"]], "MGYA1")
    expect_equal(res[["ena_accessions"]], '["ERR1","SAM1"]')
    expect_equal(res[["run_accession"]], "ERR1")
    expect_equal(res[["sample_sample_title"]], "title")
    expect_true(is.na(res[["assembly"]]))
    expect_equal(res[["publications"]], '[{"pubmed_id":1},{"pubmed_id":2}]')
    expect_equal(res[["empty_nested"]], '[{"a":1}]')
    # JSON keeps the values exactly, e.g., null and one-element arrays
    pubs <- list(list(title = "a;b", doi = NULL, authors = list("X")))
    res <- .v2_flatten(list(publications = pubs))
    expect_equal(
        res[["publications"]], '[{"title":"a;b","doi":null,"authors":["X"]}]')
    expect_identical(jsonlite::parse_json(res[["publications"]]), pubs)
    expect_false(any(startsWith(names(res), "downloads")))
    # Records with different fields are combined, values are not changed,
    # and types are converted only if the values stay the same
    df <- .v2_items_to_df(list(
        list(accession = "A", count = 1L, id = "0012", x = ""),
        list(accession = "B", count = 2L, id = "12", x = "not provided"),
        list(downloads = list(list(url = "a")))), prefix = "s")
    expect_equal(dim(df), c(3L, 4L))
    expect_equal(colnames(df), c("s_accession", "s_count", "s_id", "s_x"))
    expect_type(df[["s_count"]], "integer")
    expect_equal(df[["s_id"]], c("0012", "12", NA))
    # Empty values are NA, but other values are kept
    expect_equal(df[["s_x"]], c(NA, "not provided", NA))
    # Columns without values are dropped, except the required ones
    df <- data.frame(
        analysis_accession = "MGYA1", run_accession = NA, a = NA, b = "x")
    expect_equal(
        colnames(.v2_drop_empty_metadata(df)),
        c("analysis_accession", "run_accession", "b"))
    expect_equal(.v2_convert_lossless(c("1.10", "2")), c("1.10", "2"))
    expect_equal(.v2_convert_lossless(c("0.5", NA)), c(0.5, NA))
    expect_equal(.v2_convert_lossless(c("TRUE", "FALSE")), c(TRUE, FALSE))
    expect_equal(.v2_convert_lossless(c("7", "NA")), c("7", "NA"))
    expect_equal(nrow(.v2_items_to_df(list())), 0L)
    expect_equal(dim(.v2_items_to_df(list(list(), list()))), c(2L, 0L))
    # Analyses have always the same columns
    df <- .v2_analyses_to_df(list(list(
        accession = "MGYA1", study_accession = "MGYS1",
        pipeline_version = "V5", run = NULL)))
    expect_equal(df[["analysis_accession"]], "MGYA1")
    expect_true(is.na(df[["run_accession"]]))
    expect_true(all(c("assembly_accession", "instrument_model") %in%
        colnames(df)))
    # Metadata keeps its path, so nothing clashes or is dropped
    df <- .v2_items_to_df(list(list(
        accession = "S1", sample_title = "a",
        metadata = list(sample_title = "b", lat = 1))))
    expect_equal(
        colnames(df),
        c("accession", "sample_title", "metadata_sample_title",
            "metadata_lat"))
    expect_equal(df[["sample_title"]], "a")
    expect_equal(df[["metadata_sample_title"]], "b")
    expect_equal(nrow(.v2_items_to_df(list(list(metadata = list())))), 1L)
})

test_that("Pagination", {
    # Mock API with 250 items, 100 per page
    local_mocked_bindings(.v2_get_json = function(urls){
        lapply(urls, function(url){
            page <- as.integer(sub(".*[?&]page=([0-9]+).*", "\\1", url))
            idx <- seq((page - 1L) * 100L + 1L, min(page * 100L, 250L))
            list(count = 250L, items = lapply(idx, function(i)
                list(accession = paste0("A", i))))
        })
    })
    expect_equal(searchStudies()[["accession"]], paste0("A", 1:250))
    expect_equal(nrow(searchStudies(max.hits = 120)), 120L)
    expect_equal(searchStudies(max.hits = 5)[["accession"]], paste0("A", 1:5))
    res <- searchSamples(c("MGYS1", "MGYS2"), max.hits = 300)
    expect_equal(as.vector(table(res[["study_accession"]])), c(250L, 50L))
    # Failed page makes the whole result fail
    local_mocked_bindings(.v2_get_json = function(urls){
        lapply(urls, function(url){
            if( grepl("page=2", url) ) return(NULL)
            list(count = 150L, items = list(list(accession = "A")))
        })
    })
    expect_message(res <- searchStudies(), "could not be fetched completely")
    expect_equal(nrow(res), 0L)
})

test_that("Pagination uses the page size given by the API", {
    # Mock API that gives at most 40 items per page, although 100 are asked
    local_mocked_bindings(.v2_get_json = function(urls){
        lapply(urls, function(url){
            page <- as.integer(sub(".*[?&]page=([0-9]+).*", "\\1", url))
            idx <- seq((page - 1L) * 40L + 1L, min(page * 40L, 250L))
            list(count = 250L, items = lapply(idx, function(i)
                list(accession = paste0("A", i))))
        })
    })
    expect_equal(searchStudies()[["accession"]], paste0("A", 1:250))
    expect_equal(
        searchStudies(max.hits = 120)[["accession"]], paste0("A", 1:120))
    # A first page without items falls back to the asked page size
    local_mocked_bindings(.v2_get_json = function(urls){
        lapply(urls, function(url) list(count = 0L, items = list()))
    })
    expect_equal(nrow(searchStudies()), 0L)
})

test_that("searchStudies orders by accession and drops duplicates", {
    # Mock API that ignores the order and lists A2 on two pages, like MGnify
    # does when the order is not given
    urls_seen <- character(0)
    local_mocked_bindings(.v2_get_json = function(urls){
        urls_seen <<- c(urls_seen, urls)
        lapply(urls, function(url){
            page <- as.integer(sub(".*[?&]page=([0-9]+).*", "\\1", url))
            ids <- list(c("A1", "A2"), c("A2", "A3"))[[page]]
            list(count = 4L, items = lapply(ids, function(id)
                list(accession = id)))
        })
    })
    expect_warning(res <- searchStudies(), "1 duplicated")
    expect_equal(res[["accession"]], c("A1", "A2", "A3"))
    expect_true(all(grepl("order=accession", urls_seen, fixed = TRUE)))
    urls_seen <- character(0)
    suppressWarnings(searchStudies(order = "-updated_at"))
    expect_true(all(grepl("order=-updated_at", urls_seen, fixed = TRUE)))
})

test_that("Mixed types and failed studies", {
    local_mocked_bindings(.v2_get_all = function(paths, ...){
        list(
            list(list(accession = "S1", metadata = list(depth = 5))),
            list(list(accession = "S2", metadata = list(depth = "5 m"))),
            NULL)
    })
    expect_message(
        res <- searchSamples(c("MGYS1", "MGYS2", "MGYS3")), "MGYS3")
    expect_equal(res[["metadata_depth"]], c("5", "5 m"))
    expect_equal(res[["study_accession"]], c("MGYS1", "MGYS2"))
})

# Uses MGnify. Skipped if MGnify is not available.
test_online("searchStudies, searchSamples and searchAnalyses", {
    res <- searchStudies(biome = "root:Engineered", max.hits = 5)
    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), 5L)
    expect_true(all(startsWith(res[["biome_lineage"]], "root:Engineered")))
    res <- searchStudies(search = "MGYS00005116")
    expect_true("MGYS00005116" %in% res[["accession"]])
    res <- searchStudies(order = "-accession", max.hits = 3)
    expect_equal(
        res[["accession"]], sort(res[["accession"]], decreasing = TRUE))

    res <- searchSamples("MGYS00005116")
    expect_true(all(res[["study_accession"]] == "MGYS00005116"))
    expect_true(all(c("SAMN02594912", "SAMN02595037") %in% res[["accession"]]))

    res <- searchAnalyses(
        "MGYS00005116", samples = TRUE, biosamples = TRUE, studies = TRUE,
        analyses = TRUE, assemblies = TRUE, qc = TRUE)
    expect_true(nrow(res) >= 6L)
    expect_true(all(c(
        "analysis_accession", "pipeline_version", "study_biome_lineage",
        "analysis_quality_control_summary_after_filtering_total_reads",
        "analysis_submitted_count", "analysis_qc_statistics_sequence_count",
        "sample_sample_title", "biosamples_isolation_source") %in%
        colnames(res)))
    # By default, only the fields of the analysis lists are there
    res_default <- searchAnalyses("MGYS00005116")
    expect_equal(sort(res_default[["analysis_accession"]]),
        sort(res[["analysis_accession"]]))
    expect_true("sample_sample_title" %in% colnames(res_default))
    expect_false(any(grepl("^(analysis|study|biosamples)_", setdiff(
        colnames(res_default), c("analysis_accession", "study_accession")))))
    # The other sources do not change them
    cols <- colnames(res_default)
    expect_true(all(cols %in% colnames(res)))
    rows <- match(res_default[["analysis_accession"]],
        res[["analysis_accession"]])
    expect_equal(res_default[, cols], res[rows, cols], ignore_attr = TRUE)
    expect_true(all(c("V5", "V6") %in% res[["pipeline_version"]]))
    expect_equal(
        unique(res[["study_biome_lineage"]]),
        "root:Engineered:Built environment")
    # Quality control of pipeline v4 from its files
    v4 <- res[["pipeline_version"]] == "V4.1"
    expect_true(all(res[v4, "analysis_qc_statistics_sequence_count"] > 0))
    # Study-specific attributes come from BioSamples
    res <- searchAnalyses("MGYS00001502", biosamples = TRUE)
    expect_true(all(c("biosamples_donor_status", "biosamples_genotype") %in%
        colnames(res)))
    # Missing studies give a message, not an error
    expect_message(
        res <- searchAnalyses(c("MGYS00005116", "MGYS99999999")),
        "MGYS99999999")
    expect_true(nrow(res) > 0L)
})

test_that("Read quality control files", {
    path <- tempfile()
    writeLines(c(
        "Submitted nucleotide sequences\t16492",
        "Nucleotide sequences after undetermined bases filtering\t16144"),
        path)
    expect_equal(
        .v2_read_qc_file(path),
        c("Submitted nucleotide sequences" = "16492",
            "Nucleotide sequences after undetermined bases filtering" =
                "16144"))
    count <- tempfile()
    writeLines("112110", count)
    expect_equal(.v2_read_qc_file(count), c("112110") |> setNames(""))
    expect_equal(.v2_read_qc_file(NA_character_), character(0))
    # Only the files of each pipeline are read. Analyses without a results
    # folder have none.
    details <- list(
        list(results_dir = "https://x/ERR1_FASTQ"),
        list(results_dir = "https://x/ERR2_FASTQ"),
        list())
    urls_seen <- character(0)
    local_mocked_bindings(.v2_download = function(urls, ...){
        urls_seen <<- urls
        ifelse(grepl("count$", urls), count, path)
    })
    res <- .v2_qc_files(details, c("V4.1", "V5", "V5"))
    expect_equal(sort(urls_seen), sort(c(
        "https://x/ERR1_FASTQ/ERR1_FASTQ.fasta.submitted.count",
        "https://x/ERR1_FASTQ/qc-statistics/summary.out",
        "https://x/ERR2_FASTQ/qc_summary",
        "https://x/ERR2_FASTQ/qc-statistics/summary.out")))
    expect_equal(nrow(res), 3L)
    expect_equal(res[["analysis_submitted_count"]], c(112110L, NA, NA))
    expect_equal(
        res[["analysis_qc_summary_Submitted nucleotide sequences"]],
        c(NA, 16492L, NA))
    expect_equal(dim(.v2_qc_files(list(list()), "V5")), c(1L, 0L))
})

test_that("Sample attributes from BioSamples", {
    sample <- list(
        accession = "SAM1", name = "x",
        characteristics = list(
            `donor status` = list(list(text = "PD")),
            age = list(list(
                text = "7", unit = "weeks", ontologyTerms = list("UO_1"))),
            `body site` = list(list(text = "feces"), list(text = "colon")),
            name = list(list(text = "y")),
            bmi = list(list(text = "Missing: Not provided"))),
        `_links` = list(self = list(href = "a")))
    record <- .v2_biosample_record(sample)
    expect_false("accession" %in% names(record))
    expect_equal(record[["donor status"]], "PD")
    expect_equal(record[["age"]], "7")
    expect_equal(record[["age_unit"]], "weeks")
    expect_equal(record[["age_ontologyTerms"]], "UO_1")
    expect_equal(record[["body site"]], "feces;colon")
    expect_equal(record[["name"]], "x")
    expect_equal(record[["characteristics_name"]], "y")
    df <- .v2_items_to_df(list(record), prefix = "biosamples")
    expect_equal(df[["biosamples_bmi"]], "Missing: Not provided")
    expect_equal(df[["biosamples_age"]], 7L)
    expect_false(any(grepl("_links", colnames(df))))
    expect_length(.v2_biosample_record(NULL), 0L)
})

test_that("Downloaded files are read from the cache", {
    cache <- tempfile()
    old <- options(MGnifyR.cache = cache)
    on.exit(options(old))
    downloaded <- character(0)
    local_mocked_bindings(.v2_multi_download = function(urls, files){
        downloaded <<- c(downloaded, urls)
        # Write the URL into the file, so that the content can be checked
        for( i in seq_along(urls) ) writeLines(urls[[i]], files[[i]])
        data.frame(success = TRUE, status_code = rep(200, length(urls)))
    })
    urls <- c("https://x/a.json", "https://x/b.json", "https://x/a.json")
    # First time, each URL is downloaded once, also if given twice
    paths <- .v2_download(urls, ext = ".json")
    expect_equal(downloaded, c("https://x/a.json", "https://x/b.json"))
    expect_equal(paths[[1L]], paths[[3L]])
    expect_equal(unname(vapply(paths, readLines, character(1))), urls)
    # Second time, the files come from the cache, and only new URLs are
    # downloaded. Each file type is cached the same way.
    downloaded <- character(0)
    paths2 <- .v2_download(urls, ext = ".json")
    expect_length(downloaded, 0L)
    expect_equal(paths2, paths)
    for( ext in c(".tsv", ".txt") ){
        .v2_download("https://x/c", ext = ext)
        downloaded <- character(0)
        path <- .v2_download("https://x/c", ext = ext)
        expect_length(downloaded, 0L)
        expect_true(endsWith(path, ext))
        expect_equal(readLines(path), "https://x/c")
    }
    # The same URL with a different extension is a different file
    expect_false(.v2_download("https://x/c", ".tsv") ==
        .v2_download("https://x/c", ".txt"))
})

test_that("Robust downloads and messages", {
    # The login token is sent only to MGnify API
    seen <- list()
    local_mocked_bindings(multi_download = function(urls, files, ...){
        seen[[length(seen) + 1L]] <<- list(urls = urls, args = list(...))
        data.frame(success = TRUE, status_code = 200)
    })
    .v2_state[["token"]] <- "secret"
    res <- .v2_multi_download(
        c(paste0(.V2_URL, "studies/"), paste0(.BIOSAMPLES_BULK_URL, "SAM1")),
        c("a", "b"))
    rm("token", envir = .v2_state)
    expect_equal(res[["status_code"]], c(200, 200))
    api <- Filter(function(x) startsWith(x$urls, .V2_URL), seen)[[1L]]
    other <- Filter(
        function(x) startsWith(x$urls, .BIOSAMPLES_BULK_URL), seen)[[1L]]
    expect_equal(api$args$httpheader, "Authorization: Bearer secret")
    expect_null(other$args$httpheader)
    # Downloads have timeouts, so that a stalled connection does not hang
    expect_equal(api$args$connecttimeout, .V2_CONNECT_TIMEOUT)
    expect_equal(other$args$low_speed_time, .V2_STALL_TIMEOUT)

    # A broken JSON file gives NULL and is removed from the cache
    broken <- tempfile(fileext = ".json")
    writeLines("<html>error</html>", broken)
    local_mocked_bindings(.v2_download = function(urls, ...) broken)
    expect_equal(.v2_get_json("https://x"), list(NULL))
    expect_false(file.exists(broken))

    # Messages show at most 10 accessions
    expect_equal(.v2_shorten(c("A", "B")), "'A', 'B'")
    expect_equal(.v2_shorten(LETTERS, 2L), "'A', 'B', ...")
    expect_message(
        .message_not_found(list(NULL, 1, NULL), c("A", "B", "C"), "Studies"),
        "Studies could not be fetched \\(2\\): 'A', 'C'.*try again later")

    # BioSamples attributes do not overwrite other fields
    sample <- list(name = "x", characteristics = list(
        age = list(list(text = "7", unit = "weeks")),
        age_unit = list(list(text = "days"))))
    record <- .v2_biosample_record(sample)
    expect_equal(record[["age_unit"]], "weeks")
    expect_equal(record[["characteristics_age_unit"]], "days")
    expect_equal(record[["name"]], "x")
})

test_that("samples are fetched from the sample lists of studies", {
    # Samples that are not in the lists are fetched by themselves
    local_mocked_bindings(
        .v2_get_all = function(paths, ...) list(list(
            list(accession = "S1", sample_title = "a"),
            list(accession = "S2", sample_title = "b"))),
        .v2_get_json = function(urls){
            expect_equal(urls, paste0(.V2_URL, "samples/S3"))
            list(list(accession = "S3", sample_title = "c", studies = list()))
        })
    res <- .v2_get_samples(c("S2", "S3", "S1"), "MGYS1")
    expect_equal(vapply(res, `[[`, "", "sample_title"), c("b", "c", "a"))
    # The details of a sample have the same fields as the lists
    expect_false("studies" %in% names(res[[2L]]))
})

test_that("BioSamples are fetched in batches of 100", {
    ids <- paste0("SAMEA", seq_len(150))
    local_mocked_bindings(.v2_get_json = function(urls){
        expect_length(urls, 2L)
        expect_true(all(startsWith(urls, .BIOSAMPLES_BULK_URL)))
        batches <- sub(.BIOSAMPLES_BULK_URL, "", urls, fixed = TRUE) |>
            strsplit(",", fixed = TRUE)
        expect_equal(lengths(batches), c(100L, 50L))
        # The results are named, to check that the names do not end up in
        # the names of the samples. SAMEA2 is not found.
        res <- lapply(batches, function(batch){
            batch <- setdiff(batch, "SAMEA2")
            setNames(lapply(batch, function(id) list(accession = id)), batch)
        })
        setNames(res, urls)
    })
    res <- .v2_get_biosamples(ids)
    expect_length(res, 150L)
    expect_null(res[[2L]])
    expect_equal(res[[150L]][["accession"]], "SAMEA150")
    expect_null(names(res))
})


test_that("Each source of information is fetched only on request", {
    called <- character(0)
    add <- function(name){
        function(analyses, ...){
            called <<- c(called, name)
            analyses
        }
    }
    local_mocked_bindings(
        .v2_list_analyses = function(study.ids) data.frame(
            analysis_accession = "MGYA1", study_accession = "MGYS1",
            sample_accession = "SAM1"),
        .v2_add_analysis_info = add("analyses"),
        .v2_add_assembly_info = add("assemblies"),
        .v2_add_study_info = add("studies"),
        .v2_add_sample_info = add("samples"),
        .v2_add_biosample_info = add("biosamples"))
    # By default, only the analysis lists are fetched
    searchAnalyses("MGYS1")
    expect_length(called, 0L)
    searchAnalyses("MGYS1", biosamples = TRUE, assemblies = TRUE)
    expect_setequal(called, c("biosamples", "assemblies"))
    # Quality control files need the details of the analyses
    called <- character(0)
    searchAnalyses("MGYS1", qc = TRUE)
    expect_equal(called, "analyses")
})

test_that("Failed downloads are tried again, except missing files", {
    attempts <- list()
    local_mocked_bindings(.v2_multi_download = function(urls, files){
        attempts[[length(attempts) + 1L]] <<- urls
        # "a" fails once, "b" is not found, and "c" works
        n <- length(attempts)
        ok <- urls == "c" | (urls == "a" & n > 1L)
        for( f in files[ok] ) writeLines("x", f)
        status <- ifelse(urls == "b", 404, ifelse(ok, 200, 500))
        data.frame(success = TRUE, status_code = status)
    })
    files <- tempfile(c("a", "b", "c"))
    .v2_download_retry(c("a", "b", "c"), files)
    # The second attempt has only the failed file, and the missing file is
    # not tried again
    expect_equal(attempts, list(c("a", "b", "c"), "a"))
    expect_equal(file.exists(files), c(TRUE, FALSE, TRUE))
    # Broken downloads are not left in the cache
    expect_false(any(file.exists(paste0(files, ".part"))))
})

test_that("Details of assemblies", {
    analyses <- data.frame(
        analysis_accession = c("MGYA1", "MGYA2", "MGYA3"),
        study_accession = "MGYS1",
        assembly_accession = c("ERZ1", NA, "ERZ1"))
    local_mocked_bindings(.v2_get_json = function(urls){
        # Each assembly is fetched once
        expect_equal(urls, .v2_url("assemblies/ERZ1"))
        list(list(
            accession = "ERZ1", run_accession = "ERR1;ERR2",
            assembler_name = "metaspades"))
    })
    res <- .v2_add_assembly_info(analyses)
    expect_equal(res[["assembly_accession"]], c("ERZ1", NA, "ERZ1"))
    expect_equal(
        res[["assembly_run_accession"]], c("ERR1;ERR2", NA, "ERR1;ERR2"))
    expect_equal(
        res[["assembly_assembler_name"]], c("metaspades", NA, "metaspades"))
    # Analyses without assemblies are not changed
    expect_identical(.v2_add_assembly_info(analyses[2L, ]), analyses[2L, ])
})

test_that("Queries of searches", {
    queries <- list()
    local_mocked_bindings(.v2_get_all = function(
            paths, query = list(), max.hits = NULL){
        queries[[length(queries) + 1L]] <<- c(list(path = paths), query)
        list(list(list(accession = "A1")))
    })
    searchStudies(pipeline.version = "5.0", biome = "root:Engineered")
    expect_equal(queries[[1L]][["has_analyses_from_pipeline"]], "V5")
    expect_equal(queries[[1L]][["biome_lineage"]], "root:Engineered")
    # Samples of all studies
    res <- searchSamples(biome = "root:Host-associated", search = "gut")
    expect_equal(queries[[2L]][["path"]], "samples/")
    expect_equal(queries[[2L]][["search"]], "gut")
    expect_equal(res[["accession"]], "A1")
})

test_that("Failed quality control files are reported", {
    local_mocked_bindings(.v2_download = function(urls, ...){
        rep(NA_character_, length(urls))
    })
    expect_message(
        res <- .v2_qc_files(list(list(results_dir = "https://x/ERR1")), "V5"),
        "Could not download 2 quality control files")
    expect_equal(nrow(res), 1L)
})
