#' Fetch taxonomic and functional profiles
#'
#' @details
#' \code{fetchStudies} fetches all the analyses of the given studies, and
#' \code{fetchAnalyses} fetches specific analyses, e.g., a filtered result of
#' \code{\link{searchAnalyses}}. Columns are analyses and rows are taxa or,
#' in functional profiles, functions, e.g., GO terms. Taxa are matched across
#' studies by their taxonomy lineage.
#'
#' \code{type} tells which profile is fetched. Only analyses that have the
#' type are fetched.
#' \itemize{
#'   \item Taxonomic profiles: \code{"SSU"} and \code{"LSU"} (small and
#'   large subunit rRNA), \code{"mOTUs"} (marker genes of metagenomes and
#'   metatranscriptomes, pipelines v5-v6), \code{"UNITE"} and
#'   \code{"ITSoneDB"} (fungal ITS), \code{"PR2"} (eukaryotes),
#'   \code{"DADA2-SILVA"} and \code{"DADA2-PR2"}
#'   (amplicon sequence variants, pipeline v6), and \code{"contigs"}
#'   (number of contigs per taxon in assemblies, pipeline v6). Pipelines
#'   v1-v3 have a single profile, which is fetched with \code{"SSU"}.
#'   \item Functional profiles of metagenomes, metatranscriptomes and
#'   assemblies: \code{"InterPro"}, \code{"GO"} and \code{"GO-slim"}
#'   (pipelines v1-v6), \code{"Pfam"} and \code{"KO"} (KEGG orthologs,
#'   pipelines v5-v6), and \code{"antiSMASH"} (biosynthetic gene clusters of
#'   assemblies, pipelines v5-v6). In pipeline v6, only assemblies have
#'   functional profiles.
#' }
#'
#' Profiles of different pipeline versions are based on different reference
#' databases, so they might not be directly comparable. To fetch a single
#' pipeline version, subset the result of \code{\link{searchAnalyses}} by
#' \code{pipeline_version}. Data is cached for the R session; to keep it
#' between sessions, set \code{options(MGnifyR.cache = "path/to/cache")}.
#'
#' @param study.ids A character vector specifying study accessions, e.g.,
#' \code{"MGYS00005116"}, or a \code{data.frame} returned by
#' \code{\link{searchStudies}}.
#'
#' @param analysis.ids A character vector specifying analysis accessions, e.g.,
#' \code{"MGYA00383350"}, or a \code{data.frame} returned by
#' \code{\link{searchAnalyses}}. Giving the \code{data.frame} is faster.
#'
#' @param type A character vector specifying the profiles to fetch, e.g.,
#' \code{"SSU"} or \code{c("SSU", "GO-slim")}. See details for the types.
#' (By default: \code{type = "SSU"})
#'
#' @inheritParams searchStudies
#'
#' @details
#' The arguments \code{samples}, \code{biosamples}, \code{studies},
#' \code{analyses}, \code{assemblies} and \code{qc} choose the metadata in
#' \code{colData}, as in \code{\link{searchAnalyses}}. In
#' \code{fetchAnalyses}, they are used only if \code{analysis.ids} are
#' accessions; for a \code{data.frame}, give them to
#' \code{\link{searchAnalyses}} instead.
#'
#' @return
#' For a single type, a \code{TreeSummarizedExperiment} with
#' \itemize{
#'   \item \code{assay}: \code{"counts"}
#'   \item \code{colData}: the columns of \code{\link{searchAnalyses}}, i.e.,
#'   the metadata of the analyses, their studies and samples
#'   \item \code{rowData}: taxonomy ranks, or the descriptions of the rows of
#'   functional profiles
#' }
#'
#' For multiple types, a \code{MultiAssayExperiment} with a
#' \code{TreeSummarizedExperiment} for each type that has data. The profiles
#' are linked by sample, e.g., the SSU profile of an amplicon analysis and
#' the GO profile of a metagenomic analysis of the same sample.
#'
#' \code{NULL} is returned if no data could be fetched.
#'
#' @examples
#' # Fetch all analyses of a study
#' tse <- fetchStudies("MGYS00005116")
#' tse
#'
#' # Select analyses and fetch them
#' analyses <- searchAnalyses("MGYS00005116")
#' analyses <- analyses[analyses[["pipeline_version"]] == "V5", ]
#' tse <- fetchAnalyses(analyses)
#'
#' \donttest{
#' # Fetch a functional profile of a study
#' tse <- fetchStudies("MGYS00005610", type = "GO-slim")
#'
#' # Fetch taxonomic and functional profiles into MultiAssayExperiment
#' mae <- fetchStudies("MGYS00005610", type = c("SSU", "GO-slim"))
#'
#' # Browse studies and fetch them
#' studies <- searchStudies(biome = "root:Host-associated:Human", max.hits = 5)
#' tse <- fetchStudies(studies)
#' }
#'
#' @seealso
#' \code{\link[MGnifyR:searchStudies]{searchStudies}}
#'
#' @name fetchStudies
NULL

#' @rdname fetchStudies
#' @export
fetchStudies <- function(
        study.ids, type = "SSU", samples = FALSE, biosamples = FALSE,
        studies = FALSE, analyses = FALSE, assemblies = FALSE, qc = FALSE){
    study.ids <- .v2_ids(study.ids, "study.ids", "accession")
    .check_type(type)
    # Fetching studies is the same as fetching all their analyses
    res <- searchAnalyses(
        study.ids, samples = samples, biosamples = biosamples,
        studies = studies, analyses = analyses, assemblies = assemblies,
        qc = qc)
    no_analyses <- setdiff(study.ids, res[["study_accession"]])
    if( length(no_analyses) > 0L ){
        message(
            "No analyses found for studies: '",
            paste(no_analyses, collapse = "', '"), "'")
    }
    res <- .v2_fetch_analyses(res, type)
    return(res)
}

#' @rdname fetchStudies
#' @export
fetchAnalyses <- function(
        analysis.ids, type = "SSU", samples = FALSE, biosamples = FALSE,
        studies = FALSE, analyses = FALSE, assemblies = FALSE, qc = FALSE){
    .check_type(type)
    .check_sources(
        samples = samples, biosamples = biosamples, studies = studies,
        analyses = analyses, assemblies = assemblies, qc = qc)
    # The input is either a table from searchAnalyses() or accessions. The
    # table is used as is, which is fastest. For accessions, search the
    # analyses to get the same table, which needs one request per analysis.
    res <- analysis.ids
    if( !is.data.frame(res) ){
        ids <- .v2_ids(analysis.ids, "analysis.ids", "analysis_accession")
        res <- .v2_analyses_by_id(
            ids, samples = samples, biosamples = biosamples,
            studies = studies, analyses = analyses, assemblies = assemblies,
            qc = qc)
    }
    # Check that the table has the columns that are needed for fetching
    required <- names(.V2_ANALYSIS_COLS)
    if( nrow(res) > 0L && !all(required %in% colnames(res)) ){
        stop(
            "'analysis.ids' must be a character vector of accessions or a ",
            "data.frame returned by searchAnalyses().", call. = FALSE)
    }
    res <- .v2_fetch_analyses(res, type)
    return(res)
}

################################ HELP FUNCTIONS ################################

# Fetch the profiles of analyses. analyses is the table returned by
# searchAnalyses(), and type is checked already by the caller.
.v2_fetch_analyses <- function(analyses, type){
    if( nrow(analyses) == 0L ){
        warning("No analyses found.", call. = FALSE)
        return(NULL)
    }
    # Fetch each analysis only once
    duplicated_ids <- duplicated(analyses[["analysis_accession"]])
    analyses <- analyses[!duplicated_ids, , drop = FALSE]
    # A single type is returned as TreeSE
    if( length(type) == 1L ){
        res <- .v2_fetch_tse(analyses, type)
        return(res)
    }
    # Multiple types are fetched one at a time and combined into MAE. Types
    # without data are left out.
    type <- .v2_type_names(type) |> unique()
    tses <- lapply(type, function(x) .v2_fetch_tse(analyses, x))
    names(tses) <- type
    is_empty <- vapply(tses, is.null, logical(1))
    tses <- tses[!is_empty]
    if( length(tses) == 0L ){
        warning("No data could be fetched.", call. = FALSE)
        return(NULL)
    }
    res <- .v2_build_mae(tses, analyses)
    return(res)
}

# The types of profiles and the names of their tables. "study" is the name of
# the summary table of a study, and "analysis" is the name of the table of an
# analysis. The names are matched ignoring case. For example, SSU tables are
#   study:    "ERP1_taxonomy_abundances_SSU_v5.0.tsv" (pipelines v4-v5),
#             "ERP1_taxonomy_abundances_v3.0.tsv" (pipelines v1-v3) or
#             "ERP1_SILVA-SSU_study_summary.tsv" (pipeline v6)
#   analysis: "ERR1_MERGED_FASTQ_SSU_OTU.tsv" (pipelines v4-v5),
#             "ERR1_MERGED_FASTQ_otu.tsv" (pipelines v1-v3) or
#             "SILVA-SSU" (pipeline v6)
# Phylum-level summary tables ("ERP1_phylum_taxonomy_abundances_...") are not
# matched, because the names must start with the study accession and "_". NA
# means that there are no such tables: DADA2 tables of analyses are not
# available, and contig tables of analyses list contigs instead of counts.
# "functional" tells whether the rows are taxa or, e.g., GO terms. "ranks"
# are the taxonomy ranks of the reference database. rowData always has them,
# so that data of different studies has the same columns.
.V2_TYPES <- data.frame(
    type = c(
        "SSU", "LSU", "mOTUs", "UNITE", "ITSoneDB", "PR2", "DADA2-SILVA",
        "DADA2-PR2", "contigs",
        "InterPro", "GO", "GO-slim", "Pfam", "KO", "antiSMASH"),
    study = c(
        paste0(
            "^[^_]+_taxonomy_abundances_(SSU_)?v[0-9.]+\\.tsv$",
            "|_SILVA-SSU_study_summary\\.tsv$"),
        paste0(
            "^[^_]+_taxonomy_abundances_LSU_v[0-9.]+\\.tsv$",
            "|_SILVA-LSU_study_summary\\.tsv$"),
        "_motus_study_summary\\.tsv$",
        paste0(
            "^[^_]+_taxonomy_abundances_unite_v[0-9.]+\\.tsv$",
            "|_UNITE_study_summary\\.tsv$"),
        paste0(
            "^[^_]+_taxonomy_abundances_itsonedb_v[0-9.]+\\.tsv$",
            "|_ITSoneDB_study_summary\\.tsv$"),
        "^[^_]+_PR2_study_summary\\.tsv$",
        "_DADA2-SILVA_.*study_summary\\.tsv$",
        "_DADA2-PR2_.*study_summary\\.tsv$",
        "^[^_]+_taxonomy_study_summary\\.tsv$",
        "_IPR_abundances_|_interpro_study_summary\\.tsv$",
        "_GO_abundances_|_go_study_summary\\.tsv$",
        "_GO-slim_abundances_|_goslim_study_summary\\.tsv$",
        "_pfam_study_summary\\.tsv$",
        "_ko_study_summary\\.tsv$",
        "_antismash_study_summary\\.tsv$"),
    analysis = c(
        "_SSU_OTU\\.tsv$|(FASTQ|FASTA|SFF)_otu\\.tsv$|^SILVA-SSU$",
        "_LSU_OTU\\.tsv$|^SILVA-LSU$",
        "\\.mOTUs\\.tsv$|^motus$",
        "_ITS_UNITE_OTU\\.tsv$|^UNITE$",
        "_ITS_ITSoneDB_OTU\\.tsv$|^ITSoneDB$",
        "^PR2$",
        NA,
        NA,
        NA,
        "_InterPro\\.csv$|_interpro_summary\\.tsv",
        "_GO\\.csv$|_go_summary\\.tsv",
        "_GO_slim\\.csv$|_goslim_summary\\.tsv",
        "_pfam\\.tsv$|_pfam_summary\\.tsv",
        "_ko\\.tsv$|_ko_summary\\.tsv",
        "_antismash_summary\\.(csv|tsv)"),
    functional = rep(c(FALSE, TRUE), c(9L, 6L)),
    ranks = c(
        rep("superkingdom,kingdom,phylum,class,order,family,genus,species", 2L),
        "kingdom,phylum,class,order,family,genus,species",
        rep("superkingdom,kingdom,phylum,class,order,family,genus,species", 2L),
        paste0(
            "domain,supergroup,division,subdivision,class,order,family,genus,",
            "species"),
        "superkingdom,kingdom,phylum,class,order,family,genus,species",
        paste0(
            "domain,supergroup,division,subdivision,class,order,family,genus,",
            "species"),
        "domain,kingdom,phylum,class,order,family,genus,species",
        rep(NA_character_, 6L)))

.check_type <- function(type){
    known <- .V2_TYPES[["type"]]
    is_character <- is.character(type) && length(type) > 0L &&
        !anyNA(type)
    if( !is_character || !all(tolower(type) %in% tolower(known)) ){
        stop(
            "'type' must be one or more of: '",
            paste(known, collapse = "', '"), "'.", call. = FALSE)
    }
    return(NULL)
}

# Search analyses by analysis accessions. The API has no search by analysis
# accessions, so the studies of the analyses are found first, and the given
# analyses are taken from the analysis lists of the studies. Information
# from other sources, e.g., BioSamples, is added only for the given analyses,
# because a study can have many more analyses. ... are the arguments of
# searchAnalyses() that choose the sources.
.v2_analyses_by_id <- function(ids, ...){
    details <- paste0("analyses/", ids) |> .v2_url() |> .v2_get_json()
    .message_not_found(details, ids, "Analyses")
    study_ids <- lapply(details, `[[`, "study_accession") |> unlist() |>
        unique()
    if( length(study_ids) == 0L ){
        return(data.frame())
    }
    listed <- .v2_list_analyses(study_ids)
    # Keep the given analyses in the given order. IDs that were not found are
    # dropped.
    ids <- ids[ids %in% listed[["analysis_accession"]]]
    rows <- match(ids, listed[["analysis_accession"]])
    analyses <- listed[rows, , drop = FALSE]
    # The sample lists of the studies need one request per 100 samples of
    # the studies. If fewer samples are needed, they are fetched one by one.
    n_needed <- unique(analyses[["sample_accession"]]) |> length()
    n_pages <- vapply(study_ids, function(id){
        n_samples <- listed[listed[["study_accession"]] == id,
            "sample_accession"] |> unique() |> length()
        n <- ceiling(n_samples / 100)
        return(n)
    }, numeric(1)) |> sum()
    analyses <- .v2_add_info(
        analyses, ..., sample.lists = n_needed > n_pages)
    return(analyses)
}

# Fetch one type of profile, e.g., "SSU", of analyses. Returns a TreeSE, or
# NULL if no analysis has the type. analyses is the table of analyses returned
# by searchAnalyses().
#
# Each study has a summary table for each pipeline version. It has all the
# analyses of the study, so it is the fastest way to get the data. However,
# the summary tables are not always complete, e.g., in pipeline v6. The
# missing analyses are fetched from their own tables.
.v2_fetch_tse <- function(analyses, type){
    type <- .v2_type_names(type)
    info <- .V2_TYPES[.V2_TYPES[["type"]] == type, ]
    # Analyses are referred to by their accessions in all the steps below
    rownames(analyses) <- analyses[["analysis_accession"]]
    # Study details list the files of the study, including the summary tables
    study_ids <- unique(analyses[["study_accession"]])
    studies <- paste0("studies/", study_ids) |> .v2_url() |> .v2_get_json()
    names(studies) <- study_ids
    # Amplicons are marker genes, so they have no functional profiles. Do not
    # look for their tables, as it would be slow for large studies.
    candidates <- analyses
    if( info[["functional"]] ){
        is_amplicon <- analyses[["experiment_type"]] %in%
            c("Amplicon", "Metabarcoding")
        candidates <- analyses[!is_amplicon, , drop = FALSE]
    }

    # Step 1: Get the tables. "sources" tells the table (file name) of each
    # analysis that was found.
    res <- .v2_collect_tables(candidates, studies, info)
    tables <- res[["tables"]]
    sources <- res[["sources"]]

    # Step 2: Tell which analyses have no data, and give up if none has
    found <- rownames(analyses) %in% names(sources)
    .message_not_found_tables(rownames(analyses)[!found], type)
    if( !any(found) ){
        warning("No ", type, " data could be fetched.", call. = FALSE)
        return(NULL)
    }

    # Step 3: Add the table name to the metadata of the analyses (colData)
    analyses <- .v2_add_table_info(analyses, tables, sources)

    # Step 4: Build the TreeSE. Rows of functional profiles are described by
    # the descriptions from the tables, and taxa by their lineage.
    descriptions <- NULL
    if( info[["functional"]] ){
        descriptions <- .v2_combine_row_data(res[["row_data"]])
    }
    ranks <- strsplit(info[["ranks"]], ",", fixed = TRUE)[[1L]]
    tse <- .v2_build_tse(
        tables, analyses, row_data = descriptions, ranks = ranks)
    .message_mixed(
        tse[["pipeline_version"]], "pipeline versions",
        "Their profiles are based on different reference databases.")
    return(tse)
}

# Get the tables of analyses. Analyses are first taken from the summary tables
# of the studies, and the remaining analyses one by one from their own tables.
# Returns a list like .v2_study_tables().
.v2_collect_tables <- function(analyses, studies, info){
    res <- list(tables = list(), sources = character(0), row_data = list())
    if( nrow(analyses) > 0L ){
        res <- .v2_study_tables(analyses, studies, info)
    }
    missing <- setdiff(rownames(analyses), names(res[["sources"]]))
    if( length(missing) > 0L && !is.na(info[["analysis"]]) ){
        more <- .v2_analysis_tables(missing, info)
        for( name in names(res) ){
            res[[name]] <- c(res[[name]], more[[name]])
        }
    }
    return(res)
}

# Keep the analyses that are in the tables, in the order of the columns, and
# add the table name to their metadata
.v2_add_table_info <- function(analyses, tables, sources){
    ids <- lapply(tables, colnames) |> unlist()
    analyses <- analyses[ids, , drop = FALSE]
    analyses[["abundance_table"]] <- sources[ids]
    return(analyses)
}

# The types as they are written in .V2_TYPES, e.g., "mOTUs" for "motus"
.v2_type_names <- function(type){
    known <- .V2_TYPES[["type"]]
    rows <- match(tolower(type), tolower(known))
    res <- known[rows]
    return(res)
}

# Combine TreeSEs of different types into MAE. Columns of the TreeSEs are
# analyses. They are linked by the sample that was analysed, e.g., the SSU
# profile of an amplicon analysis and the GO profile of a metagenomic
# analysis of the same sample. A sample can have several analyses in a
# TreeSE, e.g., from different pipeline versions.
#' @importFrom MultiAssayExperiment MultiAssayExperiment ExperimentList
.v2_build_mae <- function(tses, analyses){
    # The sample of each analysis. Analyses without a sample are their own
    # samples.
    samples <- analyses[["sample_accession"]]
    no_sample <- is.na(samples)
    samples[no_sample] <- analyses[["analysis_accession"]][no_sample]
    names(samples) <- analyses[["analysis_accession"]]
    # The links between the samples and the columns of each TreeSE
    sample_map <- lapply(names(tses), function(type){
        ids <- colnames(tses[[type]])
        df <- data.frame(assay = type, primary = samples[ids], colname = ids)
        return(df)
    })
    sample_map <- do.call(rbind, sample_map)
    sample_map[["assay"]] <- factor(sample_map[["assay"]], names(tses))
    rownames(sample_map) <- NULL
    # Sample metadata, i.e., the columns that are not about the analyses.
    # The details of analyses and assemblies ("analysis_*" and "assembly_*")
    # are in colData of the TreeSEs.
    analysis_cols <- setdiff(
        names(.V2_ANALYSIS_COLS), c("study_accession", "sample_accession"))
    is_analysis_col <- colnames(analyses) %in% analysis_cols |
        grepl("^(analysis|assembly)_", colnames(analyses))
    col_data <- analyses[, !is_analysis_col, drop = FALSE]
    col_data <- col_data[!duplicated(samples), , drop = FALSE]
    rownames(col_data) <- samples[!duplicated(samples)]
    col_data <- col_data[unique(sample_map[["primary"]]), , drop = FALSE]
    col_data <- .v2_drop_empty_metadata(col_data)
    mae <- MultiAssayExperiment(
        experiments = ExperimentList(tses), colData = DataFrame(col_data),
        sampleMap = DataFrame(sample_map))
    return(mae)
}

# Tell which analyses are not in the result, because they have no table or
# their table could not be downloaded, e.g., "No GO-slim abundance table
# could be fetched for 2 analyses ...". At most 10 of the analyses are shown.
.message_not_found_tables <- function(ids, type){
    if( length(ids) == 0L ){
        return(NULL)
    }
    message(
        "No ", type, " abundance table could be fetched for ", length(ids),
        " analyses (no such table, or it could not be downloaded): ",
        .v2_shorten(ids))
    return(NULL)
}

# Give a message if the data has multiple values of something that should not
# be mixed, e.g., pipeline versions
.message_mixed <- function(x, what, why){
    x <- x[!is.na(x)] |> unique()
    if( length(x) > 1L ){
        values <- paste(x, collapse = ", ")
        message(
            "The data includes analyses from multiple ", what, " (",
            values, "). ", why, " Consider subsetting ",
            "the result of searchAnalyses() and fetching the analyses with ",
            "fetchAnalyses().")
    }
    return(NULL)
}

# Get the analyses from the summary tables of studies. info is the row of
# .V2_TYPES of the type. Returns a list with
# - "tables": abundance tables whose columns are analyses
# - "sources": the file name of the table of each analysis that was found,
#   e.g., c(MGYA1 = "S_taxonomy_abundances_SSU_v5.0.tsv")
# - "row_data": the descriptions of the rows of functional tables
.v2_study_tables <- function(analyses, studies, info){
    # Each study has its own summary tables for each pipeline version. Below,
    # such a combination of study and pipeline version is called a group.
    groups <- unique(analyses[, c("study_accession", "pipeline_version")])

    # Step 1: Find the URLs of the tables of each group from study downloads.
    # A group can have multiple tables, e.g., DADA2 tables of different
    # amplified regions.
    group_urls <- lapply(seq_len(nrow(groups)), function(i){
        study <- studies[[groups[i, "study_accession"]]]
        urls <- .v2_select_tables(
            study[["downloads"]], groups[i, "pipeline_version"], info)
        return(urls)
    })
    urls <- unlist(group_urls) |> as.character() |> .v2_fix_summary_urls()
    url_group <- rep(seq_len(nrow(groups)), lengths(group_urls))

    # Step 2: Download all the tables at once
    paths <- .v2_download_tables(urls)

    # Step 3: Read the tables, and keep the columns of the analyses of the
    # group of each table
    read <- if( info[["functional"]] ) .v2_read_functional_table else
        .v2_read_table
    tables <- list()
    sources <- character(0)
    row_data <- list()
    for( i in seq_along(urls) ){
        file <- basename(urls[[i]])
        tab <- .v2_try_read(paths[[i]], read, file)
        if( is.null(tab) ){
            next
        }
        # Functional tables also have the descriptions of the rows
        if( is.list(tab) ){
            row_data <- c(row_data, list(tab[["row_data"]]))
            tab <- tab[["counts"]]
        }
        group <- groups[url_group[[i]], ]
        tab <- .v2_group_columns(tab, analyses, group, names(sources))
        # The table might not have any of the analyses
        if( ncol(tab) == 0L ){
            next
        }
        tables <- c(tables, list(tab))
        sources[colnames(tab)] <- file
    }
    res <- list(tables = tables, sources = sources, row_data = row_data)
    return(res)
}

# Rename the columns of a summary table from runs or assemblies to the
# analyses of its group, i.e., study and pipeline version. An analysis that is
# already in another table (done) is not taken again.
.v2_group_columns <- function(tab, analyses, group, done){
    in_group <-
        analyses[["study_accession"]] %in% group[["study_accession"]] &
        analyses[["pipeline_version"]] %in% group[["pipeline_version"]]
    is_new <- !rownames(analyses) %in% done
    group_analyses <- analyses[in_group & is_new, , drop = FALSE]
    keys <- .v2_column_keys(group_analyses, tab)
    tab <- .v2_match_columns(tab, keys)
    # Some summary tables have no counts for an analysis, even though the
    # analysis has counts in its own table. Such analyses are fetched from
    # their own tables instead.
    tab <- .v2_drop_empty_columns(tab)
    return(tab)
}

# Drop the columns without counts. Summary tables of older pipelines have the
# unclassified reads in "Root". A column with counts only in "Root" has no
# taxa, so it is dropped too.
.v2_drop_empty_columns <- function(tab){
    is_taxon <- rownames(tab) != "Root"
    counts <- tab[is_taxon, , drop = FALSE] |> colSums()
    tab <- tab[, counts > 0, drop = FALSE]
    return(tab)
}

# Workaround for MGnify API: some URLs of pipeline v6 summary tables lack the
# "study-summaries" folder where the tables are. This can be removed when the
# API is fixed.
.v2_fix_summary_urls <- function(urls){
    in_folder <- grepl("/study-summaries/", urls, fixed = TRUE)
    urls[!in_folder] <- sub(
        "/([^/]+_study_summary\\.tsv)$", "/study-summaries/\\1",
        urls[!in_folder])
    return(urls)
}

# Get the tables of individual analyses. This is slower than using the
# summary tables of studies, so it is used only for the analyses that are
# missing from the summary tables. info is the row of .V2_TYPES of the type.
# Returns a list like .v2_study_tables().
.v2_analysis_tables <- function(ids, info){
    message(
        "Fetching ", length(ids), " analyses one by one, because they are ",
        "missing from the summary tables of studies. This can take a while.")
    # Step 1: Find the URL of the table of each analysis from its downloads
    details <- paste0("analyses/", ids) |> .v2_url() |> .v2_get_json()
    urls <- vapply(details, function(detail){
        url <- .v2_select_analysis_table(detail[["downloads"]], info)
        return(url)
    }, character(1))

    # Step 2: Download all the tables at once
    paths <- .v2_download_tables(urls)

    # Step 3: Read the tables. Each table has one column, the analysis.
    read <- .v2_read_analysis_table
    if( info[["functional"]] ){
        read <- .v2_read_functional_analysis_table
    }
    tables <- list()
    sources <- character(0)
    row_data <- list()
    for( i in seq_along(ids) ){
        file <- basename(urls[[i]])
        tab <- .v2_try_read(paths[[i]], read, file)
        if( is.null(tab) ){
            next
        }
        # Functional tables also have the descriptions of the rows
        if( is.list(tab) ){
            row_data <- c(row_data, list(tab[["row_data"]]))
            tab <- tab[["counts"]]
        }
        colnames(tab) <- ids[[i]]
        tables <- c(tables, list(tab))
        sources[ids[[i]]] <- file
    }
    res <- list(tables = tables, sources = sources, row_data = row_data)
    return(res)
}

# Download tables. Returns local paths, NA for tables that have no URL or that
# could not be downloaded. Tables can be listed by the API but missing from
# the FTP server, so tell which ones failed.
.v2_download_tables <- function(urls){
    has_url <- !is.na(urls)
    paths <- rep(NA_character_, length(urls))
    paths[has_url] <- .v2_download(urls[has_url], ext = ".tsv")
    failed <- urls[has_url & is.na(paths)] |> unique()
    if( length(failed) > 0L ){
        # Show at most 3 of the URLs. Some tables are listed by the API but
        # missing from the FTP server, and the server can also fail for a
        # while, so tell the user how to get the missing data in that case.
        message(
            "Could not download ", length(failed), " tables: ",
            .v2_shorten(failed, 3L), ". They might be missing from the ",
            "server, or the server did not respond. If the server did not ",
            "respond, try again later; only the missing tables are ",
            "downloaded again.")
    }
    return(paths)
}

# Read a table if it was downloaded. Returns NULL for an empty table, which
# has no data. A table that cannot be read is skipped, so that one broken
# table does not prevent fetching the other data.
.v2_try_read <- function(path, read, name = basename(path)){
    if( is.na(path) ){
        return(NULL)
    }
    tab <- tryCatch(read(path), error = function(e){
        message("Could not read ", name, ": ", conditionMessage(e))
        NULL
    })
    return(tab)
}

# Get a field, e.g., "url", of each download listed by the API
.v2_download_field <- function(downloads, name){
    res <- vapply(downloads, function(d){
        value <- d[[name]] %||% NA_character_
        return(value)
    }, character(1))
    return(res)
}

# Select the URLs of the summary tables of a type and pipeline version from
# study downloads. info is the row of .V2_TYPES of the type. The tables are
# matched by their names ("alias"). Download groups look like
# "study_summary.v5.0.taxonomic_analysis_ssu_rrna" or "study_summary.v6.go".
.v2_select_tables <- function(downloads, pipeline_version, info){
    groups <- .v2_download_field(downloads, "download_group")
    aliases <- .v2_download_field(downloads, "alias")
    # Pipeline version of each download group, e.g., "5.0"
    group_version <- sub(
        "^study_summary\\.v([0-9]+(\\.[0-9]+)*)\\..*$", "\\1", groups)
    group_version <- .v2_norm_version(group_version)
    pipeline_version <- .v2_norm_version(pipeline_version)
    is_summary <- startsWith(groups, "study_summary.")
    is_version <- group_version %in% pipeline_version
    is_type <- grepl(info[["study"]], aliases, ignore.case = TRUE)
    urls <- .v2_download_field(downloads, "url")
    urls <- urls[is_summary & is_version & is_type]
    return(urls)
}

# Select the URL of the table of a type from the downloads of an analysis.
# The table is matched by its name ("alias") instead of the URL, as the files
# can be named differently, e.g., "ERR1_MERGED_FASTQ_pfam.tsv" is
# "ERR1_MERGED_FASTQ.summary.pfam". Returns NA if there is no table.
.v2_select_analysis_table <- function(downloads, info){
    aliases <- .v2_download_field(downloads, "alias")
    is_type <- grepl(info[["analysis"]], aliases, ignore.case = TRUE)
    urls <- .v2_download_field(downloads, "url")
    url <- urls[is_type][1L]
    return(url)
}

# Get the possible column names of each analysis in a summary table. Returns
# a named list, e.g., list(MGYA1 = "ERR1", MGYA2 = c("ERR2", "ERR3")).
.v2_column_keys <- function(analyses, tab){
    # Columns are usually named by the run or assembly of the analysis. An
    # analysis can have multiple runs, which MGnify separates in different
    # ways, e.g., "ERR2,ERR3", "ERR2, ERR3" or "ERR2;ERR3".
    keys <- lapply(rownames(analyses), function(id){
        ids <- c(
            analyses[id, "assembly_accession"], analyses[id, "run_accession"])
        ids <- ids[!is.na(ids)]
        ids <- strsplit(ids, "[,;]") |> unlist() |> trimws()
        return(ids)
    })
    names(keys) <- rownames(analyses)

    # Older pipelines use other names, e.g., sample IDs "ERS095011_T" or
    # legacy assembly IDs "OGXD01". They are the start of the names of the
    # result files of the analysis, e.g., "ERS095011_T_SFF" or
    # "OGXD01_FASTA_I5.tsv.gz". Fetch the analysis details to get the file
    # names, but only for the analyses that are not found otherwise.
    found <- vapply(keys, function(key){
        is_found <- any(key %in% colnames(tab))
        return(is_found)
    }, logical(1))
    not_found <- !found
    if( !any(not_found) ){
        return(keys)
    }
    ids <- names(keys)[not_found]
    details <- paste0("analyses/", ids) |> .v2_url() |> .v2_get_json()
    for( i in seq_along(ids) ){
        download_urls <- .v2_download_field(
            details[[i]][["downloads"]], "url")
        files <- c(details[[i]][["results_dir"]], download_urls) |>
            basename()
        # E.g., "ERS095011_T" from "ERS095011_T_SFF"
        legacy <- sub("_(FASTA|FASTQ|SFF).*$", "", files) |> unique()
        keys[[ids[[i]]]] <- c(keys[[ids[[i]]]], legacy)
    }
    return(keys)
}

# Rename the columns of a summary table to analysis accessions. keys is a
# named list of possible column names of each analysis. Columns that do not
# belong to the analyses are dropped.
.v2_match_columns <- function(tab, keys){
    # Column of each analysis: the first of its names that is in the table
    cols <- vapply(keys, function(key){
        key <- key[key %in% colnames(tab)]
        col <- NA_character_
        if( length(key) > 0L ){
            col <- key[[1L]]
        }
        return(col)
    }, character(1))
    # A column is ambiguous if the table has several columns with its name,
    # or if it is the column of several analyses. E.g., in pipeline v1, the
    # assembly and the reads of a run are both named by the run. These
    # analyses are not taken from the summary table, but they are fetched
    # from their own tables.
    duplicated_names <- colnames(tab)[duplicated(colnames(tab))]
    shared_cols <- cols[duplicated(cols)]
    is_ambiguous <- cols %in% c(duplicated_names, shared_cols)
    cols[is_ambiguous] <- NA_character_
    found <- !is.na(cols)
    tab <- tab[, cols[found], drop = FALSE]
    colnames(tab) <- names(keys)[found]
    return(tab)
}
