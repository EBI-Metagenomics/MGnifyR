#' Search MGnify studies, samples and analyses
#'
#' @details
#' Browse studies with \code{searchStudies} and fetch their data with
#' \code{\link{fetchStudies}}. To select specific analyses, e.g., by sample
#' metadata or pipeline version, filter the result of \code{searchAnalyses}
#' and give it to \code{\link{fetchAnalyses}}.
#'
#' \code{searchAnalyses} returns one row per analysis with the accessions,
#' pipeline version, experiment type, sequencing instrument, and the other
#' fields of the analysis in the list of analyses of its study. The other
#' sources are added with separate arguments, as each needs its own requests.
#' By default, none of them is added, which is the fastest. The data is
#' cached, so a later search with more sources fetches only the new data:
#' \itemize{
#'   \item \code{samples}: sample metadata in MGnify (\code{sample_*}),
#'   one request per 100 samples
#'   \item \code{biosamples}: all attributes of the samples in
#'   \href{https://www.ebi.ac.uk/biosamples/}{BioSamples}
#'   (\code{biosamples_*}), including study-specific variables, e.g.,
#'   \code{biosamples_disease_status}, one request per 100 samples
#'   \item \code{studies}: details and publications of the studies
#'   (\code{study_*}), two requests per study
#'   \item \code{analyses}: details of the analyses (\code{analysis_*}),
#'   e.g., the quality control of pipeline v6, one request per analysis
#'   \item \code{assemblies}: details of the assemblies
#'   (\code{assembly_*}), one request per assembly
#'   \item \code{qc}: quality control files of pipelines v1-v5
#'   (\code{analysis_*}), several files per analysis, and the details of
#'   the analyses
#' }
#' These are the same columns as in the
#' \code{colData} of \code{\link{fetchAnalyses}}. Values are as in the
#' database. Lists, e.g., the publications of a study
#' (\code{study_publications}), are JSON text that can be read with
#' \code{jsonlite::fromJSON()}.
#'
#' In MGnify, a study has samples, a sample has runs, i.e., sequencings of
#' the sample, and an analysis is a run or an assembly processed with an
#' MGnify pipeline. A sample can have multiple runs, e.g., an amplicon and a
#' metagenome, and a run can be analysed with multiple pipeline versions, so
#' a sample can have multiple analyses. An analysis usually has one run, but
#' some have multiple runs, e.g., \code{"ERR1, ERR2"}, and analyses of
#' assemblies can have no run in \code{run_accession}.
#'
#' Results are cached for the R session; to keep them between sessions, set
#' \code{options(MGnifyR.cache = "path/to/cache")}. To access private data,
#' log in first with \code{\link{mgnifyLogin}}.
#'
#' @param biome \code{NULL} or a single character value specifying a biome
#' lineage, e.g., \code{"root:Environmental:Aquatic"}. Results also include
#' the descendant biomes. Matching is not case-sensitive. The available
#' biomes are listed at
#' \url{https://www.ebi.ac.uk/metagenomics/browse/biomes}.
#' (By default: \code{biome = NULL})
#'
#' @param pipeline.version \code{NULL} or a single character value specifying
#' the version of the MGnify pipeline: \code{"V1"}, \code{"V2"}, \code{"V3"},
#' \code{"V4"}, \code{"V4.1"}, \code{"V5"}, \code{"V6"}, \code{"V6.1"} or
#' \code{"V6.2"}. The \code{"V"} and \code{".0"} can be omitted, e.g.,
#' \code{"5.0"} or \code{"5"}. (By default: \code{pipeline.version = NULL})
#'
#' @param search \code{NULL} or a single character value used to search
#' titles and accessions. (By default: \code{search = NULL})
#'
#' @param order \code{NULL} or a single character value specifying the order
#' of the results. For studies, \code{"accession"} or \code{"updated_at"},
#' and for samples, \code{"sample_title"} or \code{"updated_at"}. Prefix
#' with \code{"-"} for descending order, e.g., \code{"-updated_at"}.
#' (By default: \code{order = NULL})
#'
#' @param max.hits \code{NULL} or a single integer value specifying the maximum
#' number of results. If \code{NULL}, all results are returned.
#' (By default: \code{max.hits = NULL})
#'
#' @param study.ids A character vector specifying study accessions, e.g.,
#' \code{"MGYS00005116"}, or a \code{data.frame} returned by
#' \code{searchStudies}. In \code{searchSamples}, it can also be
#' \code{NULL} to search all samples. \code{biome}, \code{search} and
#' \code{order} can be used only when searching all samples.
#'
#' @param samples \code{TRUE} or \code{FALSE} specifying whether to add the
#' sample metadata in MGnify. (By default: \code{samples = FALSE})
#'
#' @param biosamples \code{TRUE} or \code{FALSE} specifying whether to add
#' the attributes of the samples in BioSamples.
#' (By default: \code{biosamples = FALSE})
#'
#' @param studies \code{TRUE} or \code{FALSE} specifying whether to add the
#' details and publications of the studies.
#' (By default: \code{studies = FALSE})
#'
#' @param analyses \code{TRUE} or \code{FALSE} specifying whether to add the
#' details of the analyses, e.g., the quality control of pipeline v6.
#' (By default: \code{analyses = FALSE})
#'
#' @param assemblies \code{TRUE} or \code{FALSE} specifying whether to add
#' the details of the assemblies. (By default: \code{assemblies = FALSE})
#'
#' @param qc \code{TRUE} or \code{FALSE} specifying whether to add the
#' quality control of pipelines v1-v5, e.g., the number of sequences before
#' and after quality filtering. It is read from separate files, which makes
#' the search slower. The details of the analyses are added too.
#' (By default: \code{qc = FALSE})
#'
#' @return A \code{data.frame} where each row is a study, sample or analysis.
#'
#' @examples
#' # Search studies of aquatic environments
#' studies <- searchStudies(biome = "root:Environmental:Aquatic", max.hits = 10)
#'
#' # Get samples and analyses of a study
#' samples <- searchSamples("MGYS00005116")
#' analyses <- searchAnalyses("MGYS00005116")
#' analyses <- analyses[analyses[["pipeline_version"]] == "V6", ]
#'
#' @seealso
#' \code{\link[MGnifyR:fetchStudies]{fetchStudies}}
#'
#' @name searchStudies
NULL

#' @rdname searchStudies
#' @export
searchStudies <- function(
        biome = NULL, pipeline.version = NULL, search = NULL, order = NULL,
        max.hits = NULL){
    .check_string(biome, "biome", null.ok = TRUE)
    .check_string(pipeline.version, "pipeline.version", null.ok = TRUE)
    .check_string(search, "search", null.ok = TRUE)
    .check_string(order, "order", null.ok = TRUE)
    .check_max_hits(max.hits)
    if( !is.null(pipeline.version) ){
        pipeline.version <- paste0("V", .v2_norm_version(pipeline.version))
    }
    # Without an order, the API does not keep the same order between pages,
    # so studies would be repeated in some pages and missing from others.
    # Accessions are unique, so ordering by them keeps the pages consistent.
    order <- order %||% "accession"
    query <- list(
        biome_lineage = biome, has_analyses_from_pipeline = pipeline.version,
        search = search, order = order)
    items <- .v2_get_all("studies/", query, max.hits)
    .message_failed(items)
    res <- .v2_items_to_df(items[[1L]])
    res <- .drop_duplicated_items(res)
    return(res)
}

#' @rdname searchStudies
#' @export
searchSamples <- function(
        study.ids = NULL, biome = NULL, search = NULL, order = NULL,
        max.hits = NULL){
    .check_string(biome, "biome", null.ok = TRUE)
    .check_string(search, "search", null.ok = TRUE)
    .check_string(order, "order", null.ok = TRUE)
    .check_max_hits(max.hits)
    # Either search samples of specific studies or all samples
    if( is.null(study.ids) ){
        query <- list(biome_lineage = biome, search = search, order = order)
        items <- .v2_get_all("samples/", query, max.hits)
        .message_failed(items)
        res <- .v2_items_to_df(items[[1L]])
        return(res)
    }
    study.ids <- .v2_ids(study.ids, "study.ids", "accession")
    if( !(is.null(biome) && is.null(search) && is.null(order)) ){
        stop(
            "'biome', 'search' and 'order' cannot be used together with ",
            "'study.ids'.", call. = FALSE)
    }
    items <- paste0("studies/", study.ids, "/samples/") |> .v2_get_all()
    .message_not_found(items, study.ids, "Studies")
    # Add the study to each sample, so that samples can be linked to studies
    for( i in seq_along(items) ){
        items[[i]] <- lapply(items[[i]], function(sample){
            c(sample, list(study_accession = study.ids[[i]]))
        })
    }
    # Combine the samples of all studies into one table
    res <- unlist(items, recursive = FALSE) |> .v2_items_to_df()
    # Keep at most max.hits samples
    res <- head(res, max.hits %||% nrow(res))
    return(res)
}

#' @rdname searchStudies
#' @export
searchAnalyses <- function(
        study.ids, samples = FALSE, biosamples = FALSE, studies = FALSE,
        analyses = FALSE, assemblies = FALSE, qc = FALSE){
    study.ids <- .v2_ids(study.ids, "study.ids", "accession")
    .check_sources(
        samples = samples, biosamples = biosamples, studies = studies,
        analyses = analyses, assemblies = assemblies, qc = qc)
    res <- .v2_list_analyses(study.ids)
    res <- .v2_add_info(
        res, samples = samples, biosamples = biosamples, studies = studies,
        analyses = analyses, assemblies = assemblies, qc = qc)
    return(res)
}

################################ HELP FUNCTIONS ################################

# Get the analyses of studies from the analysis lists of the studies, one
# request per 100 analyses. Returns a table with all the fields of the lists,
# e.g., the title of the sample (see .v2_analyses_to_df()).
.v2_list_analyses <- function(study.ids){
    items <- paste0("studies/", study.ids, "/analyses/") |> .v2_get_all()
    .message_not_found(items, study.ids, "Studies")
    res <- unlist(items, recursive = FALSE) |> .v2_analyses_to_df()
    return(res)
}

# Add information from other sources to a table of analyses, df. Each source
# needs its own requests, so the user chooses which ones are fetched (see
# searchAnalyses()). The columns are named by their source, e.g., "study_".
# The quality control files are found from the details of the analyses, so
# the details are also added with qc = TRUE. sample.lists tells whether the
# samples are taken from the sample lists of the studies, which is faster
# unless only a few samples of large studies are needed (see
# .v2_get_samples()).
.v2_add_info <- function(
        df, samples = FALSE, biosamples = FALSE, studies = FALSE,
        analyses = FALSE, assemblies = FALSE, qc = FALSE,
        sample.lists = TRUE){
    res <- df
    if( analyses || qc ){
        res <- .v2_add_analysis_info(res, qc)
    }
    if( assemblies ){
        res <- .v2_add_assembly_info(res)
    }
    if( studies ){
        res <- .v2_add_study_info(res)
    }
    if( samples ){
        res <- .v2_add_sample_info(res, sample.lists)
    }
    if( biosamples ){
        res <- .v2_add_biosample_info(res)
    }
    res <- .v2_drop_empty_metadata(res)
    return(res)
}

# Add the details of analyses ("analysis_*"), e.g., their quality control in
# pipeline v6, and if qc is TRUE, the quality control files of pipelines
# v1-v5. The files are on the FTP server, which can be slow, so they are not
# read by default.
.v2_add_analysis_info <- function(analyses, qc = FALSE){
    ids <- analyses[["analysis_accession"]]
    if( length(ids) == 0L ){
        return(analyses)
    }
    details <- paste0("analyses/", ids) |> .v2_url() |> .v2_get_json()
    .message_not_found(details, ids, "Details of analyses")
    # The accession is already in "analysis_accession"
    records <- lapply(details, function(detail){
        record <- detail[names(detail) != "accession"]
        return(record)
    })
    info <- .v2_items_to_df(records, prefix = "analysis")
    if( qc ){
        qc_info <- .v2_qc_files(details, analyses[["pipeline_version"]])
        info <- cbind(info, qc_info)
    }
    res <- .v2_add_columns(analyses, info)
    return(res)
}

# Add the details of the assemblies of analyses ("assembly_*"), e.g., the run
# that was assembled
.v2_add_assembly_info <- function(analyses){
    assembly_ids <- analyses[["assembly_accession"]]
    ids <- unique(assembly_ids[!is.na(assembly_ids)])
    if( length(ids) == 0L ){
        return(analyses)
    }
    assemblies <- paste0("assemblies/", ids) |> .v2_url() |> .v2_get_json()
    .message_not_found(assemblies, ids, "Details of assemblies")
    # The accession is already in "assembly_accession"
    records <- lapply(assemblies, function(assembly){
        record <- assembly[names(assembly) != "accession"]
        return(record)
    })
    info <- .v2_items_to_df(records, prefix = "assembly")
    # Take the assembly of each analysis, i.e., one row per analysis
    rows <- match(assembly_ids, ids)
    info <- info[rows, , drop = FALSE]
    res <- .v2_add_columns(analyses, info)
    return(res)
}

# The quality control files of pipelines v1-v5 in the results folder of each
# analysis, e.g., ".../ERR1_MERGED_FASTQ". Pipeline v6 has the quality
# control in the details of the analysis instead. "{name}" is the name of the
# folder, e.g., "ERR1_MERGED_FASTQ". The files have lines of names and values
# separated by a tab, or a single value:
# - v1-v3: "ERR1_MERGED_FASTQ_summary", e.g.,
#     Submitted nucleotide sequences <TAB> 400965
#     Nucleotide sequences after trimming and quality filtering <TAB> 382583
# - v4: "ERR1_MERGED_FASTQ.fasta.submitted.count", the number of submitted
#   sequences, e.g., "112110", and "qc-statistics/summary.out", e.g.,
#     sequence_count <TAB> 98567
# - v5: "qc_summary", e.g.,
#     Submitted nucleotide sequences <TAB> 16492
#   and "qc-statistics/summary.out"
# "label" names the columns, e.g.,
# "analysis_summary_Submitted nucleotide sequences".
.V2_QC_FILES <- data.frame(
    pipeline = c("1", "2", "3", "4", "4", "5", "5"),
    file = c(
        "{name}_summary", "{name}_summary", "{name}_summary",
        "{name}.fasta.submitted.count", "qc-statistics/summary.out",
        "qc_summary", "qc-statistics/summary.out"),
    label = c(
        "summary", "summary", "summary", "submitted_count", "qc_statistics",
        "qc_summary", "qc_statistics"))

# Read the quality control files of analyses (see .V2_QC_FILES). Each line
# becomes a column. Returns a data.frame with one row per analysis.
.v2_qc_files <- function(details, pipeline_version){
    dirs <- vapply(details, function(detail){
        dir <- detail[["results_dir"]] %||% NA_character_
        return(dir)
    }, character(1))
    # Major version of the pipeline, e.g., "4" for "V4.1"
    major <- .v2_norm_version(pipeline_version) |>
        sub(pattern = "\\..*$", replacement = "")
    # The files to read: one row per analysis and file. Only the files of
    # each pipeline are requested, because a request for a file that does
    # not exist takes several seconds.
    analysis_files <- data.frame(analysis = seq_along(dirs), pipeline = major)
    files <- merge(analysis_files, .V2_QC_FILES)
    files <- files[!is.na(dirs[files[["analysis"]]]), , drop = FALSE]
    # URLs of the files, e.g., ".../ERR1_MERGED_FASTQ/ERR1_MERGED_FASTQ_summary"
    urls <- vapply(seq_len(nrow(files)), function(i){
        dir <- dirs[[files[i, "analysis"]]]
        file <- sub("{name}", basename(dir), files[i, "file"], fixed = TRUE)
        url <- paste0(dir, "/", file)
        return(url)
    }, character(1))
    paths <- .v2_download(urls, ext = ".txt")
    failed <- urls[is.na(paths)]
    if( length(failed) > 0L ){
        message(
            "Could not download ", length(failed), " quality control files: ",
            .v2_shorten(failed), ". If the server did not respond, try ",
            "again later; only the missing files are downloaded again.")
    }
    # Collect the values of each analysis, e.g.,
    # list("summary_Submitted nucleotide sequences" = "400965", ...)
    items <- rep(list(list()), length(dirs))
    for( i in seq_along(urls) ){
        values <- .v2_read_qc_file(paths[[i]])
        if( length(values) == 0L ){
            next
        }
        label <- files[i, "label"]
        names(values) <- ifelse(
            nzchar(names(values)), paste0(label, "_", names(values)), label)
        a <- files[i, "analysis"]
        items[[a]] <- c(items[[a]], as.list(values))
    }
    res <- .v2_items_to_df(items, prefix = "analysis")
    return(res)
}

# Read a quality control file. Returns its values named by the text before
# the tab on each line, e.g., c("sequence_count" = "98567"). A value without
# a name, e.g., "112110", is named "". The values are text, as in the file.
.v2_read_qc_file <- function(path){
    if( is.na(path) ){
        return(character(0))
    }
    lines <- readLines(path, warn = FALSE)
    lines <- lines[nzchar(lines)]
    has_name <- grepl("\t", lines, fixed = TRUE)
    keys <- ifelse(has_name, sub("\t.*$", "", lines), "")
    values <- ifelse(has_name, sub("^[^\t]*\t", "", lines), lines)
    res <- setNames(values, keys)
    return(res)
}

# Add the details of the studies of analyses ("study_*"), e.g., the study
# biome, which tells the environment of the samples also when the samples
# have no biome of their own, and the publications of the studies.
.v2_add_study_info <- function(analyses){
    if( nrow(analyses) == 0L ){
        return(analyses)
    }
    study_ids <- unique(analyses[["study_accession"]])
    studies <- paste0("studies/", study_ids) |> .v2_url() |> .v2_get_json()
    .message_not_found(studies, study_ids, "Details of studies")
    publications <- paste0("studies/", study_ids, "/publications/") |>
        .v2_get_all()
    .message_not_found(publications, study_ids, "Publications of studies")
    records <- lapply(seq_along(study_ids), function(i){
        # The accession is already in "study_accession"
        record <- studies[[i]]
        record <- record[names(record) != "accession"]
        record[["publications"]] <- publications[[i]]
        return(record)
    })
    info <- .v2_items_to_df(records, prefix = "study")
    # Take the study of each analysis, i.e., one row per analysis
    rows <- match(analyses[["study_accession"]], study_ids)
    info <- info[rows, , drop = FALSE]
    res <- .v2_add_columns(analyses, info)
    return(res)
}

# Add the details of the samples of analyses from MGnify ("sample_*"). MGnify
# has a fixed set of metadata fields from ENA, e.g., "lat" and
# "collection_date". They are fetched from the sample lists of the studies,
# which needs one request per 100 samples, or if lists is FALSE, one by one
# (see .v2_get_samples()).
.v2_add_sample_info <- function(analyses, lists = TRUE){
    sample_ids <- analyses[["sample_accession"]]
    ids <- unique(sample_ids[!is.na(sample_ids)])
    if( length(ids) == 0L ){
        return(analyses)
    }
    study_ids <- unique(analyses[["study_accession"]])
    samples <- .v2_get_samples(ids, study_ids, lists)
    .message_not_found(samples, ids, "Details of samples")
    # The accession is already in "sample_accession"
    records <- lapply(samples, function(sample){
        record <- sample[names(sample) != "accession"]
        return(record)
    })
    info <- .v2_items_to_df(records, prefix = "sample")
    # Take the sample of each analysis, i.e., one row per analysis
    rows <- match(sample_ids, ids)
    info <- info[rows, , drop = FALSE]
    res <- .v2_add_columns(analyses, info)
    return(res)
}

# Add the attributes of the samples of analyses from BioSamples
# ("biosamples_*"). BioSamples has all the attributes of the samples,
# including study-specific variables, e.g., "disease status", "treatment" or
# "timepoint", which are not in MGnify.
.v2_add_biosample_info <- function(analyses){
    sample_ids <- analyses[["sample_accession"]]
    ids <- unique(sample_ids[!is.na(sample_ids)])
    if( length(ids) == 0L ){
        return(analyses)
    }
    biosamples <- .v2_get_biosamples(ids)
    .message_not_found(biosamples, ids, "Samples in BioSamples")
    records <- lapply(biosamples, .v2_biosample_record)
    info <- .v2_items_to_df(records, prefix = "biosamples")
    # Take the sample of each analysis, i.e., one row per analysis
    rows <- match(sample_ids, ids)
    info <- info[rows, , drop = FALSE]
    res <- .v2_add_columns(analyses, info)
    return(res)
}

# Get the details of samples from MGnify. The sample lists of the studies have
# the same fields as the details of each sample ("samples/<id>"), except the
# studies of the sample, which are already known. A list page has 100
# samples, so this needs about 100 times fewer requests than fetching each
# sample. A sample that is not in the lists, e.g., if a page could not be
# fetched, is fetched by itself. If lists is FALSE, all the samples are
# fetched by themselves, which needs fewer requests when only a few samples
# of large studies are needed. Returns a list with the sample of each id,
# NULL if it was not found.
.v2_get_samples <- function(ids, study_ids, lists = TRUE){
    items <- list()
    if( lists ){
        pages <- paste0("studies/", study_ids, "/samples/") |> .v2_get_all()
        items <- unlist(pages, recursive = FALSE)
    }
    accessions <- vapply(items, function(item){
        accession <- item[["accession"]] %||% NA_character_
        return(accession)
    }, character(1))
    # The first one if a sample is in many lists. Missing samples are NULL.
    samples <- items[match(ids, accessions)]
    is_missing <- vapply(samples, is.null, logical(1))
    if( any(is_missing) ){
        details <- paste0("samples/", ids[is_missing]) |> .v2_url() |>
            .v2_get_json()
        # Same fields as in the lists
        samples[is_missing] <- lapply(details, function(sample){
            if( is.null(sample) ){
                return(NULL)
            }
            sample <- sample[names(sample) != "studies"]
            return(sample)
        })
    }
    return(samples)
}

# Get samples from BioSamples. The bulk endpoint returns up to 100 samples per
# request, so this needs about 100 times fewer requests than fetching each
# sample. More than 100 samples per request gives an error. The samples have
# the same values as the samples fetched one by one, with two differences:
# - The attributes are as they were submitted, e.g., "disease_status". When
#   samples are fetched one by one, BioSamples curates them, e.g., renames
#   "disease_status" to "disease status" and removes some attributes. All
#   samples are fetched in bulk, so the names are the same for all data.
# - The attributes have no ontology terms, e.g., the EFO term of a disease.
# Returns a list with the sample of each id, NULL if it was not found.
.v2_get_biosamples <- function(ids){
    batches <- split(ids, ceiling(seq_along(ids) / 100L))
    urls <- vapply(batches, function(batch){
        url <- paste0(.BIOSAMPLES_BULK_URL, paste(batch, collapse = ","))
        return(url)
    }, character(1), USE.NAMES = FALSE)
    # Each result is an object of samples named by their accessions, e.g.,
    # {"SAMEA1": {...}, "SAMEA2": {...}}. Samples that were not found are
    # not in it. unname() makes sure that the samples are named only by
    # their accessions after unlist().
    results <- .v2_get_json(urls) |> unname()
    found <- unlist(results, recursive = FALSE)
    samples <- found[ids]
    names(samples) <- NULL
    return(samples)
}

# Get the fields of a sample from BioSamples. The attributes
# ("characteristics") have one or more values, each with a text and
# optionally a unit, ontology terms and a tag, e.g.,
#   "age": [{"text": "7", "unit": "weeks", "ontologyTerms": [...]}]
# Each attribute becomes fields named by the attribute, e.g., "age",
# "age_unit" and "age_ontologyTerms". Multiple values are collapsed with ";".
# A field whose name is taken, e.g., an attribute "name" or "age_unit", is
# prefixed with "characteristics_", so that no value is overwritten.
.v2_biosample_record <- function(sample){
    # The accession is already in "sample_accession"
    record <- sample[!names(sample) %in% c("accession", "characteristics")]
    attributes <- sample[["characteristics"]]
    for( name in names(attributes) ){
        values <- attributes[[name]]
        parts <- list(text = lapply(values, `[[`, "text"))
        for( part in c("unit", "ontologyTerms", "tag") ){
            parts[[part]] <- lapply(values, `[[`, part) |> unlist()
        }
        for( part in names(parts) ){
            # Skip the parts that the attribute does not have, e.g., unit
            if( part != "text" && length(parts[[part]]) == 0L ){
                next
            }
            field <- name
            if( part != "text" ){
                field <- paste0(name, "_", part)
            }
            if( field %in% names(record) ){
                field <- paste0("characteristics_", field)
            }
            record[[field]] <- .v2_collapse(parts[[part]])
        }
    }
    return(record)
}

# Collapse values into a single character value, e.g., c("A", "B") into
# "A;B". No values give NA.
.v2_collapse <- function(x){
    x <- unlist(x)
    if( length(x) == 0L ){
        return(NA_character_)
    }
    res <- as.character(x) |> paste(collapse = ";")
    return(res)
}

.check_string <- function(x, name, null.ok = FALSE){
    is_string <- is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
    if( !(is_string || (null.ok && is.null(x))) ){
        stop(
            "'", name, "' must be a single character value",
            ifelse(null.ok, " or NULL.", "."), call. = FALSE)
    }
    return(NULL)
}

.check_flag <- function(x, name){
    if( !(isTRUE(x) || isFALSE(x)) ){
        stop("'", name, "' must be TRUE or FALSE.", call. = FALSE)
    }
    return(NULL)
}

# Check the arguments of searchAnalyses(), fetchStudies() and fetchAnalyses()
# that choose the sources of information
.check_sources <- function(
        samples = FALSE, biosamples = FALSE, studies = FALSE, analyses = FALSE,
        assemblies = FALSE, qc = FALSE){
    .check_flag(samples, "samples")
    .check_flag(biosamples, "biosamples")
    .check_flag(studies, "studies")
    .check_flag(analyses, "analyses")
    .check_flag(assemblies, "assemblies")
    .check_flag(qc, "qc")
    return(NULL)
}

.check_max_hits <- function(max.hits){
    is_count <- is.numeric(max.hits) && length(max.hits) == 1L &&
        !is.na(max.hits) && max.hits %% 1 == 0 && max.hits > 0
    if( !(is.null(max.hits) || is_count) ){
        stop(
            "'max.hits' must be NULL or a single positive integer value.",
            call. = FALSE)
    }
    return(NULL)
}

# Get IDs from a character vector or from a column of a data.frame returned
# by search functions
.v2_ids <- function(x, name, col){
    if( is.data.frame(x) ){
        x <- x[[col]]
    }
    is_valid <- is.character(x) && length(x) > 0L && !anyNA(x) &&
        all(nzchar(x))
    if( !is_valid ){
        stop(
            "'", name, "' must be a character vector of accessions or a ",
            "data.frame with column '", col, "'.", call. = FALSE)
    }
    x <- unique(x)
    return(x)
}

# Tell if search failed, so that it is not mistaken for no results
.message_failed <- function(items){
    failed <- vapply(items, is.null, logical(1))
    if( any(failed) ){
        message(
            "The search could not be fetched completely. The query might be ",
            "invalid, or the server did not respond. If the server did not ",
            "respond, try again later.")
    }
    return(NULL)
}

# Drop items that were listed in multiple pages. This happens if the order of
# the items changes between pages, e.g., if the data is updated during the
# search. Then some items might also be missing, so warn about it.
.drop_duplicated_items <- function(df){
    if( !"accession" %in% colnames(df) ){
        return(df)
    }
    is_duplicated <- duplicated(df[["accession"]])
    if( any(is_duplicated) ){
        warning(
            "MGnify returned ", sum(is_duplicated), " duplicated results, ",
            "so some results might be missing. Try again later or with ",
            "order = 'accession'.", call. = FALSE)
        df <- df[!is_duplicated, , drop = FALSE]
        rownames(df) <- NULL
    }
    return(df)
}

# Give a message about accessions that could not be fetched, e.g.,
# "Studies could not be fetched (2): 'MGYS1', 'MGYS2'. ..." The cause is not
# known here: the accession might not exist or be private, or the server did
# not respond, e.g., during a short outage. In the last case, the values are
# NA in the result, so the user is told to try again. Failed requests are not
# cached, so only they are fetched again.
.message_not_found <- function(items, ids, type){
    failed <- vapply(items, is.null, logical(1))
    not_found <- ids[failed]
    if( length(not_found) > 0L ){
        message(
            type, " could not be fetched (", length(not_found), "): ",
            .v2_shorten(not_found), ". They might not exist or be private, ",
            "or the server did not respond. If the server did not respond, ",
            "try again later; only the missing data is fetched again.")
    }
    return(NULL)
}

# Quote values for a message, showing at most n of them, e.g.,
# "'MGYA1', 'MGYA2', ..."
#' @importFrom utils head
.v2_shorten <- function(x, n = 10L){
    res <- paste0("'", head(x, n), "'", collapse = ", ")
    if( length(x) > n ){
        res <- paste0(res, ", ...")
    }
    return(res)
}
