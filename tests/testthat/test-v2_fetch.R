context("v2_fetch")
# The row of .V2_TYPES of a type, e.g., "SSU"
info <- function(type) .V2_TYPES[.V2_TYPES[["type"]] == type, ]

test_that("Input check", {
    expect_error(fetchStudies())
    expect_error(fetchStudies(1))
    expect_error(fetchStudies(NULL))
    expect_error(fetchStudies(TRUE))
    expect_error(fetchStudies("MGYS1", type = TRUE))
    expect_error(fetchStudies("MGYS1", type = c("a", NA)))
    expect_error(fetchStudies("MGYS1", type = c("SSU", "unknown")))
    expect_error(fetchStudies("MGYS1", type = character(0)))
    expect_error(fetchStudies("MGYS1", type = "unknown"), "'type' must be")
    expect_error(fetchAnalyses())
    expect_error(fetchAnalyses(1))
    expect_error(fetchAnalyses(NULL))
    expect_error(fetchAnalyses("MGYA1", type = 1))
    # Tables must have the columns of searchAnalyses()
    expect_error(
        fetchAnalyses(data.frame(analysis_accession = "MGYA1")),
        "data.frame returned by searchAnalyses")
})

test_that("Parse taxonomy lineages", {
    lineage <- c(
        "sk__Bacteria;k__;p__Bacillota",
        "Root;k__Bacteria;p__Firmicutes;c__Bacilli",
        "d__Archaea",
        "unclassified")
    res <- .v2_parse_lineage(lineage)
    expect_equal(rownames(res), lineage)
    expect_equal(
        colnames(res),
        c("domain", "superkingdom", "kingdom", "phylum", "class"))
    expect_equal(res[["superkingdom"]], c("Bacteria", NA, NA, NA))
    expect_equal(res[["kingdom"]], c(NA, "Bacteria", NA, NA))
    expect_equal(res[["phylum"]], c("Bacillota", "Firmicutes", NA, NA))
    expect_equal(res[["domain"]], c(NA, NA, "Archaea", NA))
    expect_true(all(is.na(unlist(res[4, ]))))
    # The given ranks are always included, and PR2 ranks are parsed
    pr2 <- "d__Bacteria;sg__Terrabacteria;dv__Firmicutes;sdv__Firmicutes_X"
    res <- .v2_parse_lineage(pr2, ranks = c("domain", "supergroup", "genus"))
    expect_equal(
        colnames(res),
        c("domain", "supergroup", "division", "subdivision", "genus"))
    expect_equal(res[["subdivision"]], "Firmicutes_X")
    expect_true(is.na(res[["genus"]]))
})

test_that("Select summary table", {
    download <- function(group, alias, type = "Taxonomic analysis"){
        list(
            download_group = group, alias = alias, download_type = type,
            file_type = ifelse(grepl("csv$", alias), "csv", "tsv"),
            url = paste0("https://x/", alias))
    }
    downloads <- list(
        download(
            "study_summary.v5.0.taxonomic_analysis_ssu_rrna",
            "S_phylum_taxonomy_abundances_SSU_v5.0.tsv"),
        download(
            "study_summary.v5.0.taxonomic_analysis_ssu_rrna",
            "S_taxonomy_abundances_SSU_v5.0.tsv"),
        download(
            "study_summary.v5.0.taxonomic_analysis_lsu_rrna",
            "S_taxonomy_abundances_LSU_v5.0.tsv"),
        download(
            "study_summary.v5.0.functional_analysis",
            "S_GO_abundances_v5.0.tsv", "Functional analysis"),
        download(
            "study_summary.v6.amplicon", "S_closedref_SILVA-SSU_dwcready.csv"),
        download("study_summary.v6.amplicon",
            "S_DADA2-SILVA_18S-V1-V2_study_summary.tsv"),
        download("study_summary.v6.amplicon", "S_SILVA-SSU_study_summary.tsv"),
        download("study_summary.v6.amplicon", "S_PR2_study_summary.tsv"),
        download(
            "study_summary.v6.rawreads.silva-ssu",
            "S_silva-ssu_study_summary.tsv"),
        download(
            "study_summary.v6.rawreads.motus", "S_motus_study_summary.tsv"),
        download(
            "study_summary.v2.0.taxonomic_analysis",
            "S_taxonomy_abundances_v2.0.tsv"))
    url <- function(alias) paste0("https://x/", alias)
    expect_equal(.v2_select_tables(downloads, "V5", info("SSU")),
        url("S_taxonomy_abundances_SSU_v5.0.tsv"))
    expect_equal(.v2_select_tables(downloads, "V5", info("LSU")),
        url("S_taxonomy_abundances_LSU_v5.0.tsv"))
    expect_equal(.v2_select_tables(downloads, "V5", info("GO")),
        url("S_GO_abundances_v5.0.tsv"))
    # DADA2 and Darwin Core tables are not the same as closed-reference tables
    expect_equal(.v2_select_tables(downloads, "V6", info("PR2")),
        url("S_PR2_study_summary.tsv"))
    expect_equal(.v2_select_tables(downloads, "V6", info("DADA2-SILVA")),
        url("S_DADA2-SILVA_18S-V1-V2_study_summary.tsv"))
    expect_equal(.v2_select_tables(downloads, "V6", info("mOTUs")),
        url("S_motus_study_summary.tsv"))
    # Amplicons and raw reads have their own SSU tables
    expect_equal(.v2_select_tables(downloads, "V6", info("SSU")), url(c(
        "S_SILVA-SSU_study_summary.tsv", "S_silva-ssu_study_summary.tsv")))
    # The single table of older pipelines is SSU
    expect_equal(.v2_select_tables(downloads, "V2", info("SSU")),
        url("S_taxonomy_abundances_v2.0.tsv"))
    # Types that are not there are not fetched
    expect_length(.v2_select_tables(downloads, "V2", info("mOTUs")), 0L)
    expect_length(.v2_select_tables(downloads, "V5", info("mOTUs")), 0L)
    expect_length(.v2_select_tables(downloads, "V5", info("UNITE")), 0L)
    expect_length(.v2_select_tables(downloads, "V6", info("LSU")), 0L)
    expect_length(.v2_select_tables(downloads, "V4.1", info("SSU")), 0L)
    expect_length(.v2_select_tables(downloads, NA, info("SSU")), 0L)
    expect_length(.v2_select_tables(list(), "V5", info("SSU")), 0L)
})

test_that("Read and merge abundance tables", {
    path <- tempfile(fileext = ".tsv")
    writeLines(c(
        "#SampleID\tERR1\tERR2\tERZ3",
        "sk__Bacteria\t1\t0\t5",
        "Root;sk__Bacteria;k__;p__Bacillota\t2.0\t3.0\t"), path)
    tab <- .v2_read_table(path)
    expect_equal(
        rownames(tab), c("sk__Bacteria", "sk__Bacteria;k__;p__Bacillota"))
    expect_equal(colnames(tab), c("ERR1", "ERR2", "ERZ3"))
    # Missing values are zero
    expect_equal(unname(tab["sk__Bacteria", ]), c(1, 0, 5))
    expect_equal(unname(tab[2, "ERZ3"]), 0)

    # Columns are matched to analyses by run or assembly. Multiple runs are
    # split, and other names are used only if needed.
    analyses <- data.frame(
        analysis_accession = c("MGYA1", "MGYA2", "MGYA3", "MGYA4", "MGYA5"),
        run_accession = c("ERR2", "ERR7,ERR1", NA, "ERR9", "ERR8"),
        assembly_accession = c(NA, NA, "ERZ3", NA, NA),
        row.names = c("MGYA1", "MGYA2", "MGYA3", "MGYA4", "MGYA5"))
    local_mocked_bindings(.v2_get_json = function(urls){
        expect_equal(urls, .v2_url(c("analyses/MGYA4", "analyses/MGYA5")))
        list(
            list(results_dir = "https://x/ERR9_FASTQ"),
            list(downloads = list(list(url = "https://x/OGX5_FASTA_I5.tsv"))))
    })
    # Legacy name of MGYA5
    tab_legacy <- cbind(tab, OGX5 = c(4, 0))
    keys <- .v2_column_keys(analyses, tab_legacy)
    expect_equal(keys[["MGYA2"]], c("ERR7", "ERR1"))
    expect_equal(keys[["MGYA4"]], c("ERR9", "ERR9"))
    res <- .v2_match_columns(tab_legacy, keys)
    expect_equal(colnames(res), c("MGYA1", "MGYA2", "MGYA3", "MGYA5"))
    expect_equal(unname(res[, "MGYA1"]), unname(tab[, "ERR2"]))
    expect_equal(unname(res[, "MGYA5"]), c(4, 0))

    # Ambiguous columns are not matched: the assembly and the reads of a run
    # can have the same column name, or a table can have duplicated names
    keys_shared <- list(MGYA1 = "ERR1", MGYA2 = "ERR1", MGYA3 = "ERZ3")
    res_shared <- .v2_match_columns(tab, keys_shared)
    expect_equal(colnames(res_shared), "MGYA3")
    path_dup <- tempfile(fileext = ".tsv")
    writeLines(c("#SampleID\tERR1\tERR1", "sk__Bacteria\t1\t5"), path_dup)
    tab_dup <- .v2_read_table(path_dup)
    expect_equal(colnames(tab_dup), c("ERR1", "ERR1"))
    res_dup <- .v2_match_columns(tab_dup, list(MGYA1 = "ERR1"))
    expect_equal(ncol(res_dup), 0L)

    # Tables of different studies are merged by lineage
    tab2 <- matrix(
        c(4, 7), ncol = 1L,
        dimnames = list(c("sk__Bacteria", "sk__Archaea"), "MGYA9"))
    col_data <- data.frame(
        pipeline_version = rep("V5", 5L),
        row.names = c("MGYA1", "MGYA2", "MGYA3", "MGYA5", "MGYA9"))
    tse <- .v2_build_tse(list(res, tab2), col_data)
    expect_s4_class(tse, "TreeSummarizedExperiment")
    expect_true(is.matrix(assay(tse)))
    expect_equal(dim(tse), c(3L, 5L))
    expect_equal(colnames(tse), rownames(col_data))
    expect_equal(assay(tse)["sk__Archaea", "MGYA1"], 0)
    expect_equal(assay(tse)["sk__Bacteria", "MGYA9"], 4)
    expect_equal(assay(tse)["sk__Bacteria", "MGYA3"], 5)
    expect_equal(rowData(tse)[["phylum"]], c(NA, "Bacillota", NA))
    # The same taxon written with and without spaces is one taxon, named by
    # the spelling with spaces
    tab3 <- matrix(
        c(1, 2), ncol = 1L, dimnames = list(
            c("g__CandidatusNitrososphaera", "sk__Bacteria"), "MGYA1"))
    tab4 <- matrix(
        3, dimnames = list("g__Candidatus Nitrososphaera", "MGYA2"))
    tse3 <- .v2_build_tse(
        list(tab3, tab4), data.frame(row.names = c("MGYA1", "MGYA2")))
    expect_equal(
        rownames(tse3), c("g__Candidatus Nitrososphaera", "sk__Bacteria"))
    expect_equal(unname(assay(tse3)[1L, ]), c(1, 3))
    # Large data is sparse
    tse_sparse <- .v2_build_tse(list(res, tab2), col_data, max.dense = 1)
    expect_s4_class(assay(tse_sparse), "dgCMatrix")
    expect_equal(as.matrix(assay(tse_sparse)), assay(tse))
    # Empty table, i.e., nothing detected, gives zero counts
    empty <- matrix(numeric(0), 0L, 1L, dimnames = list(NULL, "MGYA1"))
    tse_empty <- .v2_build_tse(
        list(empty, tab2), col_data[c("MGYA1", "MGYA9"), , drop = FALSE])
    expect_equal(unname(assay(tse_empty)[, "MGYA1"]), c(0, 0))
    tse_empty <- .v2_build_tse(list(empty), col_data["MGYA1", , drop = FALSE])
    expect_equal(dim(tse_empty), c(0L, 1L))
})

test_that("Fetch data with mocked API", {
    path <- tempfile(fileext = ".tsv")
    writeLines(c(
        "taxonomy\tERR1\tERR2",
        "sk__Bacteria\t1\t2",
        "sk__Bacteria;k__;p__Bacillota\t3\t0"), path)
    study <- list(title = "Study 1", downloads = list(list(
        download_group = "study_summary.v6.amplicon",
        alias = "S_SILVA-SSU_study_summary.tsv",
        download_type = "Taxonomic analysis", file_type = "tsv",
        url = "https://x/S_SILVA-SSU_study_summary.tsv")))
    analysis <- function(acc, run, sample, pipeline = "V6"){
        list(
            accession = acc, study_accession = "MGYS1",
            pipeline_version = pipeline, experiment_type = "Amplicon",
            run = list(accession = run), sample = list(accession = sample))
    }
    local_mocked_bindings(
        # Study details list the summary table, and analyses have no tables
        .v2_get_json = function(urls) lapply(urls, function(u){
            if( grepl("/studies/", u) ) return(study)
            list(study_accession = "MGYS1", downloads = list())
        }),
        .v2_get_all = function(paths, ...){
            if( all(grepl("/analyses/", paths)) ){
                return(list(list(
                    analysis("MGYA1", "ERR1", "SAM1"),
                    analysis("MGYA2", "ERR2", "SAM2"),
                    analysis("MGYA3", "ERR1", "SAM1", "V5"))))
            }
            # Samples of the study
            if( all(grepl("/samples/", paths)) ){
                return(list(list(
                    list(accession = "SAM1", metadata = list(lat = 1)),
                    list(accession = "SAM2", metadata = list(lat = 2)))))
            }
            # Publications of the study
            list(list(list(pubmed_id = 1L)))
        },
        .v2_download = function(urls, ...) rep(path, length(urls)))
    expect_message(
        tse <- fetchStudies("MGYS1", samples = TRUE, studies = TRUE),
        "No SSU abundance table could be fetched for 1 analyses.*'MGYA3'")
    expect_s4_class(tse, "TreeSummarizedExperiment")
    expect_equal(colnames(tse), c("MGYA1", "MGYA2"))
    expect_equal(unname(assay(tse)[, "MGYA2"]), c(2, 0))
    expect_true(all(c(
        "analysis_accession", "study_accession", "pipeline_version",
        "run_accession", "abundance_table", "study_title",
        "sample_metadata_lat", "study_publications") %in%
        colnames(colData(tse))))
    expect_equal(tse[["sample_metadata_lat"]], c(1, 2))
    expect_equal(unique(tse[["study_title"]]), "Study 1")
    # Search result as input
    expect_message(tse <- fetchStudies(data.frame(accession = "MGYS1")))
    expect_equal(colnames(tse), c("MGYA1", "MGYA2"))
    # Table of analyses as input
    analyses <- searchAnalyses("MGYS1")
    expect_equal(colnames(fetchAnalyses(analyses[2, ])), "MGYA2")
    expect_warning(expect_message(res <- fetchAnalyses(analyses[3, ])))
    expect_null(res)
    # Accessions as input
    expect_equal(colnames(fetchAnalyses("MGYA2")), "MGYA2")
    # Missing taxonomy
    expect_warning(expect_message(
        res <- fetchStudies("MGYS1", type = "LSU")))
    expect_null(res)
})

# Uses MGnify. Skipped if MGnify is not available.
test_online("fetchStudies and fetchAnalyses", {
    # Study has analyses from V4.1, V5 and V6 pipelines
    analyses <- searchAnalyses("MGYS00005116")
    tse <- fetchAnalyses(analyses[analyses[["pipeline_version"]] == "V5", ])
    expect_s4_class(tse, "TreeSummarizedExperiment")
    expect_equal(ncol(tse), 2L)
    expect_equal(unique(tse[["study_accession"]]), "MGYS00005116")
    expect_equal(unique(tse[["pipeline_version"]]), "V5")
    expect_equal(tse[["analysis_accession"]], colnames(tse))
    expect_true(all(colSums(assay(tse)) > 0))
    expect_true(all(c("phylum", "genus") %in% colnames(rowData(tse))))

    # Only the given type is fetched. MGYA00571587 has no LSU table.
    expect_message(
        tse <- fetchAnalyses(
            analyses[analyses[["pipeline_version"]] == "V5", ],
            type = "LSU"),
        "No LSU abundance table could be fetched for 1 analyses")
    expect_equal(
        unique(tse[["abundance_table"]]), "taxonomy_abundances_LSU_v5.0.tsv")

    # MGYS00005116 (amplicons) has no mOTUs, but MGYS00010446 has. The
    # studies are searched together, because their metadata can have
    # different columns.
    both <- searchAnalyses(c("MGYS00005116", "MGYS00010446"))
    both <- both[both[["pipeline_version"]] == "V5" |
        both[["study_accession"]] == "MGYS00010446", ]
    expect_message(
        tse <- fetchAnalyses(both, type = "mOTUs"), "No mOTUs abundance table")
    expect_equal(unique(tse[["study_accession"]]), "MGYS00010446")
    tse <- fetchAnalyses(both, type = "SSU")
    expect_equal(ncol(tse), nrow(both))

    # Multiple studies are merged. Missing study does not stop fetching.
    expect_message(
        tse <- fetchStudies(
            c("MGYS00005116", "MGYS00006771", "MGYS99999999"),
            type = "SSU"),
        "MGYS99999999")
    expect_s4_class(tse, "TreeSummarizedExperiment")
    expect_equal(
        sort(unique(tse[["study_accession"]])),
        c("MGYS00005116", "MGYS00006771"))
    # Assemblies of pipeline v6 have contig taxonomy instead of SSU
    tse <- fetchStudies(c("MGYS00005116", "MGYS00010497"), type = "contigs")
    expect_equal(unique(tse[["study_accession"]]), "MGYS00010497")

    # Analysis IDs and search results give the same data
    tse <- fetchAnalyses(c("MGYA00383350", "MGYA01026512"))
    expect_equal(colnames(tse), c("MGYA00383350", "MGYA01026512"))
    analyses <- searchAnalyses("MGYS00005116")
    tse2 <- fetchAnalyses(analyses[analyses[["analysis_accession"]] %in%
        c("MGYA00383350", "MGYA01026512"), ])
    expect_equal(assay(tse2)[rownames(tse), colnames(tse)], assay(tse))
})

test_that("Abundance tables of individual analyses", {
    # Table of an analysis lists OTUs, and counts are summed by lineage
    path <- tempfile(fileext = ".tsv")
    writeLines(c(
        "# Constructed from biom file",
        "# OTU ID\tSSU\ttaxonomy\ttaxid",
        "1\t2.0\tsk__Bacteria\t2",
        "2\t3.0\tsk__Bacteria;k__;p__Bacillota\t1239",
        "3\t1.0\tsk__Bacteria\t2"), path)
    tab <- .v2_read_analysis_table(path)
    expect_equal(
        rownames(tab), c("sk__Bacteria", "sk__Bacteria;k__;p__Bacillota"))
    expect_equal(unname(tab[, 1L]), c(3, 3))
    # Empty table has no data
    writeLines("# OTU ID\tSSU\ttaxonomy", path)
    expect_silent(expect_null(.v2_try_read(path, .v2_read_analysis_table)))
    # Pipeline v6 tables have one column per rank, and they can be gzipped
    path <- tempfile(fileext = ".tsv")
    con <- gzfile(path, "w")
    writeLines(c(
        "# Count\tSuperkingdom\tKingdom\tPhylum",
        "88\tsk__Bacteria",
        "3.0\tsk__Bacteria\tk__\tp__Bacillota",
        "2\tsk__Bacteria"), con)
    close(con)
    tab <- .v2_read_analysis_table(path)
    expect_equal(
        rownames(tab), c("sk__Bacteria", "sk__Bacteria;k__;p__Bacillota"))
    expect_equal(unname(tab[, 1L]), c(90, 3))
    expect_null(.v2_try_read(NA_character_, .v2_read_analysis_table))
    # Pipeline v5 mOTUs tables separate ranks with "|"
    writeLines(c(
        "#mOTU\tconsensus_taxonomy\tcount",
        "ref_mOTU_1\tk__Bacteria|p__Bacillota\t4",
        "ref_mOTU_2\tk__Bacteria|p__Bacillota\t1",
        "Unmapped\tUnmapped\t14"), path)
    tab <- .v2_read_analysis_table(path)
    expect_equal(rownames(tab), c("k__Bacteria;p__Bacillota", "Unmapped"))
    expect_equal(unname(tab[, 1L]), c(5, 14))
    # Taxonomy assignments of OTUs in pipeline v1 have no counts
    writeLines(c(
        "6723 READ1\tRoot;k__Bacteria\t0.890",
        "892 READ2\tRoot;k__Bacteria\t1.000"), path)
    expect_message(
        expect_null(.v2_try_read(path, .v2_read_analysis_table)), "counts")

    download <- function(alias, file_type = "tsv"){
        list(
            alias = alias, download_type = "Taxonomic analysis",
            file_type = file_type, url = paste0("https://x/", alias))
    }
    downloads <- list(
        download("ERR1_SILVA-SSU.mseq"),
        download("ERR1_SILVA-SSU.html", "html"),
        download("DADA2-SILVA"), download("SILVA-SSU"), download("PR2"))
    expect_equal(
        .v2_select_analysis_table(downloads, info("SSU")),
        "https://x/SILVA-SSU")
    expect_equal(
        .v2_select_analysis_table(downloads, info("PR2")), "https://x/PR2")
    expect_true(is.na(.v2_select_analysis_table(downloads, info("LSU"))))
    # Older pipelines have a single table, which is SSU
    downloads <- list(
        download("ERR1_FASTQ_otu.tsv"),
        download("ERR1_FASTQ_pruned.tree", "tree"))
    expect_equal(
        .v2_select_analysis_table(downloads, info("SSU")),
        "https://x/ERR1_FASTQ_otu.tsv")
    # Lineages are written in the same way in all tables
    expect_equal(
        .v2_clean_lineage(c(
            "Root;k__Bacteria;p__Firmicutes", "k__Bacteria; p__Firmicutes",
            " Root", "sk__Bacteria;k__;p__Bacillota")),
        c("k__Bacteria;p__Firmicutes", "k__Bacteria;p__Firmicutes", "Root",
            "sk__Bacteria;k__;p__Bacillota"))
    downloads <- list(
        download("ERR1_FASTQ_SSU_MAPSeq.mseq.gz"),
        download("ERR1_FASTQ_SSU_OTU.tsv"), download("ERR1_FASTQ_LSU_OTU.tsv"),
        download("ERR1_FASTQ.mOTUs.tsv"), download("ERR1_motus_motus.html"))
    expect_equal(
        .v2_select_analysis_table(downloads, info("LSU")),
        "https://x/ERR1_FASTQ_LSU_OTU.tsv")
    expect_equal(
        .v2_select_analysis_table(downloads, info("SSU")),
        "https://x/ERR1_FASTQ_SSU_OTU.tsv")
    expect_equal(
        .v2_select_analysis_table(downloads, info("mOTUs")),
        "https://x/ERR1_FASTQ.mOTUs.tsv")
})

test_that("Pipeline v6 summary table URLs are fixed", {
    downloads <- list(list(
        download_group = "study_summary.v6.amplicon",
        alias = "S_SILVA-SSU_study_summary.tsv",
        download_type = "Taxonomic analysis", file_type = "tsv",
        url = "https://x/S/S_SILVA-SSU_study_summary.tsv"))
    analyses <- data.frame(
        analysis_accession = "MGYA1", study_accession = "MGYS1",
        pipeline_version = "V6", run_accession = "ERR1",
        assembly_accession = NA, row.names = "MGYA1")
    local_mocked_bindings(.v2_download = function(urls, ...){
        expect_equal(
            urls, "https://x/S/study-summaries/S_SILVA-SSU_study_summary.tsv")
        NA_character_
    })
    res <- .v2_study_tables(
        analyses, list(MGYS1 = list(downloads = downloads)), info("SSU"))
    expect_length(res[["tables"]], 0L)
})

test_that("Summary table without the analyses is skipped", {
    downloads <- list(list(
        download_group = "study_summary.v5.0.taxonomic_analysis_ssu_rrna",
        alias = "S_taxonomy_abundances_SSU_v5.0.tsv",
        download_type = "Taxonomic analysis", file_type = "tsv",
        url = "https://x/S_taxonomy_abundances_SSU_v5.0.tsv"))
    analyses <- data.frame(
        analysis_accession = "MGYA1", study_accession = "MGYS1",
        pipeline_version = "V5", run_accession = NA,
        assembly_accession = "ERZ1", row.names = "MGYA1")
    path <- tempfile(fileext = ".tsv")
    writeLines(c("#SampleID\tERZ2", "sk__Bacteria\t3"), path)
    local_mocked_bindings(
        .v2_download = function(urls, ...) path,
        .v2_get_json = function(urls) list(list(downloads = list())))
    res <- .v2_study_tables(
        analyses, list(MGYS1 = list(downloads = downloads)), info("SSU"))
    expect_length(res[["tables"]], 0L)
    expect_length(res[["sources"]], 0L)
})

test_that("Failed table downloads are reported", {
    local_mocked_bindings(.v2_download = function(urls, ...){
        ifelse(urls == "https://x/a.tsv", "a.tsv", NA_character_)
    })
    expect_message(
        paths <- .v2_download_tables(c("https://x/a.tsv", NA, "https://x/b")),
        "Could not download 1 tables: 'https://x/b'.*try again later")
    expect_equal(paths, c("a.tsv", NA, NA))
    expect_silent(.v2_download_tables(NA_character_))
})

test_that("Columns without counts are dropped from summary tables", {
    # In MGYS00001303, most columns of the summary table have only "Root",
    # although the analyses have counts in their own tables
    tab <- matrix(
        c(1, 0, 0, 1, 5, 3, 0, 0, 0), nrow = 3,
        dimnames = list(c("Root", "k__A", "k__B"), c("a1", "a2", "a3")))
    res <- .v2_drop_empty_columns(tab)
    expect_identical(colnames(res), "a2")
    expect_identical(rownames(res), rownames(tab))
})

test_that("Empty summary table has no data", {
    path <- tempfile(fileext = ".tsv")
    writeLines("#SampleID\tERR1\tERR2", path)
    expect_silent(expect_null(.v2_try_read(path, .v2_read_table)))
    writeLines(character(0), path)
    expect_silent(expect_null(.v2_try_read(path, .v2_read_table)))
})

test_that("Select functional tables", {
    download <- function(group, alias, url = alias){
        list(
            download_group = group, alias = alias,
            download_type = "Functional analysis", file_type = "tsv",
            url = paste0("https://x/", url))
    }
    downloads <- list(
        download("study_summary.v5.0.functional_analysis",
            "S_GO_abundances_v5.0.tsv", "GO_abundances_v5.0.tsv"),
        download("study_summary.v5.0.functional_analysis",
            "S_GO-slim_abundances_v5.0.tsv", "GO-slim_abundances_v5.0.tsv"),
        download("study_summary.v5.0.functional_analysis",
            "S_IPR_abundances_v5.0.tsv", "IPR_abundances_v5.0.tsv"),
        download("study_summary.v6.assembly.goslim",
            "S_goslim_study_summary.tsv"),
        download("study_summary.v6.assembly.ko", "S_ko_study_summary.tsv"))
    expect_equal(
        .v2_select_tables(downloads, "V5", info("GO-slim")),
        "https://x/GO-slim_abundances_v5.0.tsv")
    expect_equal(
        .v2_select_tables(downloads, "V5", info("GO")),
        "https://x/GO_abundances_v5.0.tsv")
    expect_equal(
        .v2_select_tables(downloads, "V6", info("GO-slim")),
        "https://x/S_goslim_study_summary.tsv")
    expect_equal(
        .v2_select_tables(downloads, "V6", info("KO")),
        "https://x/S_ko_study_summary.tsv")
    # KO has no summary table in pipeline v5
    expect_length(.v2_select_tables(downloads, "V5", info("KO")), 0L)
    # Tables of analyses are selected by their names, not by their URLs
    downloads <- list(
        download("all", "ERR1_MERGED_FASTQ_GO.csv", "ERR1.summary.go"),
        download("all", "ERR1_MERGED_FASTQ_GO_slim.csv", "ERR1.summary.slim"),
        download("all", "ERR1_MERGED_FASTQ_ko.tsv", "ERR1.summary.ko"),
        download("all", "ERZ1_goslim_summary.tsv.gz"))
    expect_equal(
        .v2_select_analysis_table(downloads, info("GO")),
        "https://x/ERR1.summary.go")
    expect_equal(
        .v2_select_analysis_table(downloads, info("GO-slim")),
        "https://x/ERR1.summary.slim")
    expect_equal(
        .v2_select_analysis_table(downloads, info("KO")),
        "https://x/ERR1.summary.ko")
    expect_equal(
        .v2_select_analysis_table(downloads[4], info("GO-slim")),
        "https://x/ERZ1_goslim_summary.tsv.gz")
    expect_true(
        is.na(.v2_select_analysis_table(downloads, info("Pfam"))))
})

test_that("Read functional tables", {
    # Summary table of a study
    path <- tempfile(fileext = ".tsv")
    writeLines(c(
        "GO\tdescription\tcategory\tERR1\tERR2",
        paste0(
            "GO:0000015\tphosphopyruvate hydratase complex\t",
            "cellular component\t1\t0"),
        "GO:0000150\trecombinase activity\tmolecular function\t2\t5"), path)
    res <- .v2_read_functional_table(path)
    expect_equal(
        res[["counts"]],
        matrix(c(1, 2, 0, 5), nrow = 2, dimnames = list(
            c("GO:0000015", "GO:0000150"), c("ERR1", "ERR2"))))
    expect_equal(res[["row_data"]][["description"]], c(
        "phosphopyruvate hydratase complex", "recombinase activity"))
    expect_equal(
        res[["row_data"]][["category"]],
        c("cellular_component", "molecular_function"))

    # Table of an analysis, pipeline v6
    writeLines(c(
        "ko\tdescription\tcount",
        "K07729\tputative transcriptional regulator\t68",
        "K22299\tHTH-type transcriptional regulator\t62"), path)
    res <- .v2_read_functional_analysis_table(path)
    expect_equal(
        res[["counts"]],
        matrix(c(68, 62), dimnames = list(c("K07729", "K22299"), NULL)))
    expect_true(all(is.na(res[["row_data"]][["category"]])))

    # Pipelines v3-v5. Descriptions can have commas and quotes.
    writeLines(c(
        '"8","K03446","MFS transporter, DHA2 family"',
        '"1","K10673","streptomycin 3"-kinase [EC:2.7.1.87]"'), path)
    res <- .v2_read_functional_analysis_table(path)
    expect_equal(rownames(res[["counts"]]), c("K03446", "K10673"))
    expect_equal(unname(res[["counts"]][, 1L]), c(8, 1))
    expect_equal(res[["row_data"]][["description"]], c(
        "MFS transporter, DHA2 family",
        'streptomycin 3"-kinase [EC:2.7.1.87]'))
    writeLines(c(
        '"GO:0009847","spore germination","biological_process","55"',
        '"GO:0055085","transmembrane transport","biological_process","44"'),
        path)
    res <- .v2_read_functional_analysis_table(path)
    expect_equal(unname(res[["counts"]][, 1L]), c(55, 44))
    expect_equal(
        unique(res[["row_data"]][["category"]]), "biological_process")
    # Empty table has no data
    writeLines(character(0), path)
    expect_null(.v2_read_functional_analysis_table(path))
    writeLines("ko\tdescription\tcount", path)
    expect_null(.v2_read_functional_analysis_table(path))

    # Descriptions of many tables are combined
    row_data <- list(
        data.frame(description = c("a", "b"), category = NA,
            row.names = c("IPR1", "IPR2")),
        data.frame(description = c("B", "c"), category = NA,
            row.names = c("IPR2", "IPR3")))
    combined <- .v2_combine_row_data(row_data)
    expect_equal(rownames(combined), c("IPR1", "IPR2", "IPR3"))
    expect_equal(combined[["description"]], c("a", "b", "c"))
    expect_false("category" %in% colnames(combined))
})

test_that("Fetch functional profiles with mocked API", {
    taxonomy_path <- tempfile(fileext = ".tsv")
    writeLines(c(
        "taxonomy\tERR1\tERR2\tERR3",
        "sk__Bacteria\t1\t2\t3"), taxonomy_path)
    goslim_path <- tempfile(fileext = ".tsv")
    writeLines(c(
        "GO\tdescription\tcategory\tERR1\tERR2",
        paste0(
            "GO:0000015\tphosphopyruvate hydratase complex\t",
            "cellular component\t1\t4"),
        "GO:0000150\trecombinase activity\tmolecular function\t2\t0"),
        goslim_path)
    download <- function(type, alias){
        list(
            download_group = "study_summary.v5.0.x", alias = alias,
            download_type = type, file_type = "tsv",
            url = paste0("https://x/", alias))
    }
    study <- list(title = "Study 1", downloads = list(
        download("Taxonomic analysis", "S_taxonomy_abundances_SSU_v5.0.tsv"),
        download("Functional analysis", "S_GO-slim_abundances_v5.0.tsv")))
    analysis <- function(acc, run, sample, type){
        list(
            accession = acc, study_accession = "MGYS1",
            pipeline_version = "V5", experiment_type = type,
            run = list(accession = run), sample = list(accession = sample))
    }
    requested <- character(0)
    local_mocked_bindings(
        .v2_get_json = function(urls) lapply(urls, function(u){
            requested <<- c(requested, u)
            if( grepl("/studies/", u) ) return(study)
            list(downloads = list())
        }),
        .v2_get_all = function(paths, ...){
            if( all(grepl("/analyses/", paths)) ){
                return(list(list(
                    analysis("MGYA1", "ERR1", "SAM1", "Metagenomic"),
                    analysis("MGYA2", "ERR2", "SAM2", "Metagenomic"),
                    analysis("MGYA3", "ERR3", "SAM1", "Amplicon"))))
            }
            # Samples of the study
            if( all(grepl("/samples/", paths)) ){
                return(list(list(
                    list(accession = "SAM1", metadata = list(lat = 1)),
                    list(accession = "SAM2", metadata = list(lat = 2)))))
            }
            # Publications of the study
            list(list(list(pubmed_id = 1L)))
        },
        .v2_download = function(urls, ...){
            ifelse(grepl("GO-slim", urls), goslim_path, taxonomy_path)
        })

    # Functional profile only. The amplicon has none, and its table is not
    # looked for.
    analyses <- searchAnalyses("MGYS1")
    requested <- character(0)
    expect_message(
        tse <- fetchAnalyses(analyses, type = "GO-slim"),
        "No GO-slim abundance table could be fetched for 1 analyses.*'MGYA3'")
    expect_false(any(grepl("analyses/MGYA3", requested)))
    expect_s4_class(tse, "TreeSummarizedExperiment")
    expect_equal(colnames(tse), c("MGYA1", "MGYA2"))
    expect_equal(rownames(tse), c("GO:0000015", "GO:0000150"))
    expect_equal(unname(assay(tse)[, "MGYA2"]), c(4, 0))
    expect_equal(
        rowData(tse)[["category"]],
        c("cellular_component", "molecular_function"))
    expect_equal(
        unique(tse[["abundance_table"]]), "S_GO-slim_abundances_v5.0.tsv")
    # Type is not case-sensitive
    expect_message(tse2 <- fetchStudies("MGYS1", type = "go-slim"))
    expect_equal(assay(tse2), assay(tse))
    # Taxonomic profile of the same study includes the amplicon
    expect_silent(tse <- fetchStudies("MGYS1", type = "SSU"))
    expect_equal(colnames(tse), c("MGYA1", "MGYA2", "MGYA3"))

    # Multiple types are combined into MAE. The amplicon and the metagenome
    # of SAM1 are linked by the sample. LSU has no data, so it is left out.
    expect_warning(
        mae <- fetchStudies(
            "MGYS1", type = c("SSU", "go-slim", "LSU"), samples = TRUE),
        "No LSU data could be fetched.")
    expect_s4_class(mae, "MultiAssayExperiment")
    expect_equal(names(mae), c("SSU", "GO-slim"))
    expect_s4_class(mae[["SSU"]], "TreeSummarizedExperiment")
    expect_equal(assay(mae[["GO-slim"]]), assay(tse2))
    expect_equal(rownames(colData(mae)), c("SAM1", "SAM2"))
    expect_equal(colData(mae)[["sample_metadata_lat"]], c(1, 2))
    expect_false("analysis_accession" %in% colnames(colData(mae)))
    expect_false(any(startsWith(colnames(colData(mae)), "analysis_")))
    sample_map <- MultiAssayExperiment::sampleMap(mae)
    ssu <- sample_map[sample_map[["assay"]] == "SSU", ]
    expect_equal(ssu[["colname"]], c("MGYA1", "MGYA2", "MGYA3"))
    expect_equal(ssu[["primary"]], c("SAM1", "SAM2", "SAM1"))
    # Types without any data give NULL
    expect_warning(
        res <- fetchStudies("MGYS1", type = c("LSU", "mOTUs")),
        "No data could be fetched.")
    expect_null(res)
})

# Uses MGnify. Skipped if MGnify is not available.
test_online("Fetch functional profiles", {
    # Pipeline v3 metagenomes
    analyses <- searchAnalyses("MGYS00001442")
    analyses <- analyses[analyses[["analysis_accession"]] %in%
        c("MGYA00096976", "MGYA00096977"), ]
    tse <- fetchAnalyses(analyses, type = "GO-slim")
    expect_s4_class(tse, "TreeSummarizedExperiment")
    expect_equal(sort(colnames(tse)), c("MGYA00096976", "MGYA00096977"))
    expect_true(all(startsWith(rownames(tse), "GO:")))
    expect_true(all(c("description", "category") %in% colnames(rowData(tse))))
    expect_true(all(colSums(assay(tse)) > 0))
})

test_that("fetchAnalyses checks the arguments of searchAnalyses", {
    expect_error(fetchAnalyses("MGYA1", biosamples = NA), "'biosamples' must")
    expect_error(fetchAnalyses("MGYA1", bogus = TRUE), "unused argument")
    expect_error(fetchStudies("MGYS1", studies = "no"), "'studies' must")
})

test_that("analyses given by accession get information only for themselves", {
    listed <- data.frame(
        analysis_accession = c("MGYA1", "MGYA2", "MGYA3"),
        study_accession = "MGYS1", sample_accession = c("S1", "S2", "S3"))
    local_mocked_bindings(
        .v2_get_json = function(urls) list(list(study_accession = "MGYS1")),
        .v2_list_analyses = function(study.ids) listed,
        .v2_add_info = function(df, ..., sample.lists){
            # One of three samples is needed, and the list of the study has
            # one page, so the lists are not cheaper
            expect_false(sample.lists)
            expect_equal(list(...), list(biosamples = FALSE))
            df
        })
    res <- .v2_analyses_by_id("MGYA2", biosamples = FALSE)
    expect_equal(res[["analysis_accession"]], "MGYA2")
})

test_that("Functional IDs on several rows are summed", {
    # In InterPro tables of pipeline v1, an entry can be on several rows with
    # different descriptions, as it was renamed between analyses
    path <- tempfile(fileext = ".tsv")
    writeLines(c(
        "IPR\tdescription\tERR1\tERR2",
        "IPR000008\tC2 calcium-dependent membrane targeting\t1\t0",
        "IPR000008\tC2 domain\t0\t2",
        "IPR000014\tPAS\t3\t1"), path)
    tab <- .v2_read_functional_table(path)
    expect_equal(rownames(tab[["row_data"]]), c("IPR000008", "IPR000014"))
    expect_equal(
        tab[["row_data"]]["IPR000008", "description"],
        "C2 calcium-dependent membrane targeting")
    tse <- .v2_build_tse(
        list(tab[["counts"]]), data.frame(row.names = c("ERR1", "ERR2")),
        row_data = tab[["row_data"]])
    expect_equal(rownames(tse), c("IPR000008", "IPR000014"))
    expect_equal(unname(assay(tse)["IPR000008", ]), c(1, 2))
    expect_equal(unname(assay(tse)["IPR000014", ]), c(3, 1))
})

test_that("Functional IDs are taken from the column before the description", {
    # Some tables of pipeline v4.1 have a column of counts before the IDs
    path <- tempfile(fileext = ".tsv")
    writeLines(c(
        "ERR1\tIPR\tdescription\tERR2",
        "0\tIPR000007\tTubby, C-terminal\t4",
        "2\tIPR000008\tC2 domain\t70"), path)
    tab <- .v2_read_functional_table(path)
    expect_equal(rownames(tab[["counts"]]), c("IPR000007", "IPR000008"))
    expect_equal(colnames(tab[["counts"]]), c("ERR1", "ERR2"))
    expect_equal(unname(tab[["counts"]][, "ERR1"]), c(0, 2))
    expect_equal(tab[["row_data"]][["description"]],
        c("Tubby, C-terminal", "C2 domain"))
})

test_that("Analyses with multiple runs are matched by each run", {
    # MGnify separates the runs in different ways
    analyses <- data.frame(
        assembly_accession = NA_character_,
        run_accession = c("ERR1,ERR2", "ERR3, ERR4", "ERR5;ERR6"),
        row.names = c("MGYA1", "MGYA2", "MGYA3"))
    tab <- matrix(1, 1, 3, dimnames = list("sk__Bacteria",
        c("ERR2", "ERR4", "ERR6")))
    keys <- .v2_column_keys(analyses, tab)
    expect_equal(keys[["MGYA2"]], c("ERR3", "ERR4"))
    res <- .v2_match_columns(tab, keys)
    expect_equal(colnames(res), c("MGYA1", "MGYA2", "MGYA3"))
})

test_that("Functional profiles of individual analyses", {
    # Pipeline v5 table of an analysis: count, ID and description
    path <- tempfile()
    writeLines(c(
        '"8","K03446","MFS transporter, DHA2 family"',
        '"1","K10673","streptomycin 3-kinase"'), path)
    # MGYA1 has a KO table, but MGYA2 has none
    local_mocked_bindings(
        .v2_get_json = function(urls) list(
            list(downloads = list(list(
                alias = "ERR1_MERGED_FASTQ_ko.tsv",
                url = "https://x/ERR1.summary.ko"))),
            list(downloads = list())),
        .v2_download = function(urls, ...) rep(path, length(urls)))
    expect_message(
        res <- .v2_analysis_tables(c("MGYA1", "MGYA2"), info("KO")),
        "one by one")
    expect_length(res[["tables"]], 1L)
    expect_equal(colnames(res[["tables"]][[1L]]), "MGYA1")
    expect_equal(unname(res[["tables"]][[1L]][, 1L]), c(8, 1))
    expect_equal(res[["sources"]], c(MGYA1 = "ERR1.summary.ko"))
    expect_equal(
        res[["row_data"]][[1L]]["K03446", "description"],
        "MFS transporter, DHA2 family")
})
