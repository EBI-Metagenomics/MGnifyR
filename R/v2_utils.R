# Internal helpers for MGnify API v2: requests, parsing and
# TreeSummarizedExperiment assembly.

################################## REQUESTS ####################################

# Build API URLs from paths, e.g., "studies/", and query parameters
#' @importFrom curl curl_escape
.v2_url <- function(path, query = list()){
    # Drop parameters that are not given
    is_null <- vapply(query, is.null, logical(1))
    query <- query[!is_null]
    if( length(query) == 0L ){
        url <- paste0(.V2_URL, path)
        return(url)
    }
    # E.g., "search=soil&page=2". Special characters, e.g., spaces, are
    # escaped.
    values <- query |> as.character() |> curl_escape()
    query <- paste0(names(query), "=", values, collapse = "&")
    url <- paste0(.V2_URL, path, "?", query)
    return(url)
}

# Fetch JSON responses. Returns a list of parsed responses, NULL for those
# that could not be fetched or read. A file that cannot be read, e.g., a
# broken response, is removed from the cache, so that it is fetched again
# next time.
#' @importFrom jsonlite read_json
.v2_get_json <- function(urls){
    paths <- .v2_download(urls, ext = ".json")
    res <- lapply(paths, function(path){
        if( is.na(path) ){
            return(NULL)
        }
        json <- tryCatch(read_json(path), error = function(e){
            unlink(path)
            NULL
        })
        return(json)
    })
    return(res)
}

# Fetch all items from list endpoints, e.g., "studies/" or
# "studies/MGYS00005116/analyses/". The API gives the items in pages of 100
# items, so all the pages are fetched. To be fast, the pages of all the paths
# are fetched in parallel. Returns a list with the items of each path, or NULL
# for a path that could not be fetched.
.v2_get_all <- function(paths, query = list(), max.hits = NULL){
    page_size <- 100L

    # Step 1: Fetch the first page of each path. It tells how many items
    # ("count") the path has, so we know how many pages there are.
    first <- rep(1L, length(paths))
    first_pages <- .v2_page_urls(paths, first, query, page_size) |>
        .v2_get_json()
    n_items <- vapply(first_pages, function(page){
        count <- page[["count"]] %||% 0
        count <- as.numeric(count)
        return(count)
    }, numeric(1))
    # Do not fetch more items than the user asked for
    n_items <- pmin(n_items, max.hits %||% Inf)
    n_pages <- .v2_n_pages(first_pages, n_items, page_size)

    # Step 2: Fetch the other pages of all paths at once. List them in a
    # table with one row per page. For example, if path 1 has 3 pages and
    # path 2 has 2 pages:
    #   path page
    #      1    2
    #      1    3
    #      2    2
    other <- data.frame(
        path = rep(seq_along(paths), n_pages - 1L),
        page = lapply(n_pages, function(n) seq_len(n)[-1L]) |> unlist() |>
            as.integer())
    other_pages <- .v2_page_urls(
        paths[other$path], other$page, query, page_size) |> .v2_get_json()

    # Step 3: Collect the items of each path from its pages
    res <- lapply(seq_along(paths), function(i){
        pages <- c(first_pages[i], other_pages[other$path == i])
        # If any page failed, return NULL rather than incomplete results
        failed <- vapply(pages, is.null, logical(1))
        if( any(failed) ){
            return(NULL)
        }
        # Each page has a list of items. Combine them into one list.
        items <- lapply(pages, `[[`, "items") |> unlist(recursive = FALSE)
        # The last page can have more items than max.hits allows
        items <- head(items, n_items[[i]])
        return(items)
    })
    return(res)
}

# URLs of pages, e.g., ".../studies/?page=2&page_size=100"
.v2_page_urls <- function(paths, pages, query, page_size){
    urls <- vapply(seq_along(paths), function(i){
        .v2_url(paths[[i]], c(query, page = pages[[i]], page_size = page_size))
    }, character(1))
    return(urls)
}

# Number of pages of each path, given the first pages and the number of items
# to fetch. The API can give fewer items per page than asked, e.g., currently
# at most 100. Take the actual page size from the first page, so that no pages
# are missed if the limit changes. A first page without items gives no
# information, so then the asked page size is used.
.v2_n_pages <- function(first_pages, n_items, page_size){
    page_items <- vapply(first_pages, function(page){
        n <- length(page[["items"]])
        return(n)
    }, integer(1))
    actual_size <- ifelse(page_items > 0L, page_items, page_size)
    # Each path has at least one page
    n_pages <- ceiling(n_items / actual_size)
    n_pages <- pmax(n_pages, 1L)
    return(n_pages)
}

################################## PARSING #####################################

# Convert counts from text to numbers. Values that are not numbers, e.g.,
# empty fields, become NA. Only the values that look like numbers are
# converted, so that as.numeric() does not warn about the others.
.v2_as_numeric <- function(x){
    is_number <- grepl(
        "^\\s*[-+]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][-+]?[0-9]+)?\\s*$", x)
    res <- rep(NA_real_, length(x))
    res[is_number] <- as.numeric(x[is_number])
    return(res)
}

# Normalise pipeline versions so that, e.g., "V5", "v5.0" and "5.0" match
.v2_norm_version <- function(x){
    x <- tolower(x)
    x <- sub("^v", "", x)
    x <- sub("\\.0$", "", x)
    return(x)
}

# Fields that are not metadata: "downloads" lists the files of studies and
# analyses, which the data is fetched from, and "_links" has links between
# the records of BioSamples
.V2_SKIP_FIELDS <- c("downloads", "_links")

# Convert a nested record into a named character vector. Each value is named
# by its path in the record, e.g., "biome_lineage" for
# {"biome": {"lineage": ...}}. Arrays are kept as JSON, e.g.,
# '["ERP1","PRJEB1"]' or '[{"pubmed_id":1,"title":"..."}]', so that no
# information is lost. The values are not changed. JSON null and empty arrays
# are NA.
#' @importFrom jsonlite toJSON
.v2_flatten <- function(x){
    fields <- setdiff(names(x), .V2_SKIP_FIELDS)
    vals <- lapply(fields, function(name){
        val <- x[[name]]
        # Missing value, e.g., null or empty array in JSON
        if( length(val) == 0L ){
            res <- setNames(NA_character_, name)
            return(res)
        }
        # Object, e.g., "biome": {"lineage": ...}. Flatten it recursively.
        if( is.list(val) && !is.null(names(val)) ){
            res <- .v2_flatten(val)
            names(res) <- paste0(name, "_", names(res))
            return(res)
        }
        # Array, e.g., ["ERP1", "PRJEB1"] or publications
        if( is.list(val) || length(val) > 1L ){
            res <- toJSON(val, auto_unbox = TRUE, null = "null", digits = NA)
            res <- setNames(as.character(res), name)
            return(res)
        }
        # Single value
        res <- setNames(as.character(val), name)
        return(res)
    })
    vals <- unlist(vals)
    return(vals)
}

# Convert a list of records into a data.frame with one row per record. The
# columns are named by the paths of the fields (see .v2_flatten()), and by the
# source, e.g., "study_", if prefix is given. Values are not changed, but a
# column is converted to numbers or logical values if every value stays
# exactly the same (see .v2_convert_lossless()).
#' @importFrom stats setNames
.v2_items_to_df <- function(items, prefix = NULL){
    # Step 1: Flatten each record into a named vector, e.g.,
    # c(accession = "MGYS1", biome_lineage = "root:Engineered")
    rows <- lapply(items, .v2_flatten)

    # Step 2: Put the records into a table. Records can have different
    # fields, so the columns are all the fields of all the records. Missing
    # fields are NA.
    cols <- lapply(rows, names) |> unlist() |> unique()
    if( length(cols) == 0L ){
        res <- data.frame(row.names = seq_along(items))
        return(res)
    }
    mat <- matrix(
        NA_character_, nrow = length(rows), ncol = length(cols),
        dimnames = list(NULL, cols))
    for( i in seq_along(rows) ){
        row <- rows[[i]]
        mat[i, names(row)] <- row
    }
    # An empty value, e.g., "" in sample metadata, is a missing value. Other
    # values, e.g., "not provided", are kept as they are, because they can
    # have meaning.
    mat[!is.na(mat) & mat == ""] <- NA_character_
    if( !is.null(prefix) ){
        colnames(mat) <- paste0(prefix, "_", colnames(mat))
    }
    df <- as.data.frame(mat, check.names = FALSE)

    # Step 3: Convert columns to numbers or logical values when no value
    # changes
    df[] <- lapply(df, .v2_convert_lossless)
    return(df)
}

# Drop metadata columns that have no values, i.e., they are NA for all rows,
# e.g., metadata fields that none of the samples has. The columns of
# .V2_ANALYSIS_COLS are always kept, as they are needed for fetching data.
.v2_drop_empty_metadata <- function(df){
    is_empty <- vapply(df, function(col) all(is.na(col)), logical(1))
    is_required <- colnames(df) %in% names(.V2_ANALYSIS_COLS)
    df <- df[, !is_empty | is_required, drop = FALSE]
    return(df)
}

# Convert text to numbers or logical values, but only if every value stays
# exactly the same, i.e., converting the result back to text gives the
# original text. Otherwise, the text is returned. For example, "12" and
# "0.5" become numbers, but "0012", "1.10", "" and "not provided" keep the
# column as text.
#' @importFrom utils type.convert
.v2_convert_lossless <- function(x){
    converted <- type.convert(x, as.is = TRUE)
    if( is.character(converted) ){
        return(x)
    }
    back <- as.character(converted)
    is_same <- identical(is.na(back), is.na(x)) &&
        all(back == x, na.rm = TRUE)
    res <- if( is_same ) converted else x
    return(res)
}

# The columns that analyses always have, so that search results can be used
# directly for fetching data. The names are the columns of the result, and the
# values are the fields of the API.
.V2_ANALYSIS_COLS <- c(
    analysis_accession = "accession",
    study_accession = "study_accession",
    pipeline_version = "pipeline_version",
    experiment_type = "experiment_type",
    sample_accession = "sample_accession",
    run_accession = "run_accession",
    assembly_accession = "assembly_accession",
    instrument_platform = "run_instrument_platform",
    instrument_model = "run_instrument_model"
)

# Convert analyses from the analysis lists of studies into a table. The first
# columns are those of .V2_ANALYSIS_COLS, and missing values are NA. They are
# followed by all the other fields of the lists, named by their path, e.g.,
# "sample_sample_title" and "run_instrument_model", so that all the
# information that was fetched is kept.
.v2_analyses_to_df <- function(items){
    df <- .v2_items_to_df(items)
    cols <- .V2_ANALYSIS_COLS
    res <- lapply(cols, function(col){
        if( !col %in% colnames(df) ){
            values <- rep(NA_character_, nrow(df))
        } else{
            values <- as.character(df[[col]])
        }
        return(values)
    })
    res <- as.data.frame(res)
    # The fields that are not in the columns above. Fields renamed above,
    # e.g., "accession", are not repeated.
    other <- df[, !colnames(df) %in% c(cols, names(cols)), drop = FALSE]
    res <- cbind(res, other)
    return(res)
}

# Add columns of information to a table of analyses. The table can already
# have some of the columns from the analysis lists, e.g.,
# "sample_sample_title". They are replaced by the columns of info, which are
# from the details of the same records and have all their fields. The
# columns of .V2_ANALYSIS_COLS are never replaced.
.v2_add_columns <- function(analyses, info){
    is_replaced <- colnames(analyses) %in% colnames(info) &
        !colnames(analyses) %in% names(.V2_ANALYSIS_COLS)
    info <- info[, !colnames(info) %in% names(.V2_ANALYSIS_COLS), drop = FALSE]
    res <- cbind(analyses[, !is_replaced, drop = FALSE], info)
    rownames(res) <- NULL
    return(res)
}

########################## TREESUMMARIZEDEXPERIMENT ############################

# Read a study summary table. The first column contains the taxonomy lineage
# and the rest are runs or assemblies.
#' @importFrom utils read.delim
.v2_read_table <- function(path){
    # Reading the counts directly as numbers is much faster than converting
    # them afterwards
    first <- readLines(path, n = 1L)
    # An empty table has no data
    if( length(first) == 0L ){
        return(NULL)
    }
    header <- strsplit(first, "\t", fixed = TRUE)
    n_counts <- length(header[[1L]]) - 1L
    col_classes <- c("character", rep("numeric", n_counts))
    tab <- read.delim(
        path, check.names = FALSE, comment.char = "", quote = "",
        colClasses = col_classes)
    if( nrow(tab) == 0L ){
        return(NULL)
    }
    mat <- tab[, -1L, drop = FALSE] |> as.matrix()
    # Subsetting makes duplicated column names unique, e.g., "ERR1.1". Keep
    # the original names so that duplicates can be recognised.
    colnames(mat) <- header[[1L]][-1L]
    mat[is.na(mat)] <- 0
    rownames(mat) <- .v2_clean_lineage(tab[[1L]])
    return(mat)
}

# Read the abundance table of a single analysis. The table lists taxa with
# their counts. The result is a one-column matrix with lineages as row names,
# like the summary tables. The tables have different formats, which are
# recognised from the header line:
# - pipelines v4-5: OTU, count, lineage and taxonomy ID, e.g.,
#     # OTU ID      SSU_rRNA  taxonomy                      taxid
#     5820          261.0     sk__Bacteria;k__;p__Bacillota  1239
# - pipeline v5 mOTUs: mOTU, lineage separated by "|" and count, e.g.,
#     #mOTU         consensus_taxonomy           count
#     ref_mOTU_1    k__Bacteria|p__Bacillota     4
# - pipeline v6: count and one column per rank. Lower ranks are missing if
#   the taxon is not classified to them, e.g.,
#     # Count  Superkingdom  Kingdom  Phylum
#     3        sk__Bacteria  k__      p__Bacillota
# Other files, e.g., OTU assignments without counts in pipeline v1, are not
# abundance tables and give an error.
.v2_read_analysis_table <- function(path){
    # The tables can be gzipped. readLines() decompresses them.
    lines <- readLines(path, warn = FALSE)
    header <- lines[startsWith(lines, "#")]
    data_lines <- lines[!startsWith(lines, "#") & nzchar(lines)]
    # An empty table has no data
    if( length(data_lines) == 0L ){
        return(NULL)
    }
    # Split the lines into fields, and get a field of each line
    rows <- strsplit(data_lines, "\t")
    field <- function(i) vapply(rows, function(row) row[i], character(1))

    if( any(startsWith(header, "# Count")) ){
        # Pipeline v6: the ranks are joined into a lineage
        counts <- field(1L)
        lineage <- vapply(rows, function(row){
            paste(row[-1L], collapse = ";")
        }, character(1))
    } else if( any(startsWith(header, "#mOTU")) ){
        # Pipeline v5 mOTUs
        counts <- field(3L)
        lineage <- field(2L)
        lineage <- gsub("|", ";", lineage, fixed = TRUE)
    } else{
        # Pipelines v4-5
        counts <- field(2L)
        lineage <- field(3L)
    }

    counts <- .v2_as_numeric(counts)
    if( is.na(counts) |> all() ){
        stop("The table has no counts.", call. = FALSE)
    }
    counts[is.na(counts)] <- 0
    # The same lineage can be on many lines, e.g., different OTUs of the same
    # taxon, so the counts are summed by lineage
    lineage <- .v2_clean_lineage(lineage)
    mat <- matrix(counts)
    mat <- rowsum(mat, lineage, reorder = FALSE)
    return(mat)
}

# Read a functional summary table of a study. The first column has the IDs,
# e.g., GO terms, followed by their description, the GO category in GO
# tables, and the counts of each run or assembly, e.g.,
#   GO          description         category            ERR1  ERR2
#   GO:0000015  phosphopyruvate...  cellular component  32    27
# Some tables of pipeline v4.1 have a column of counts before the IDs, e.g.,
#   ERR1  IPR        description        ERR2
#   0     IPR000007  Tubby, C-terminal  4
# The IDs are therefore taken from the column before the description, or
# from the first column if there is no description.
# Returns a list with "counts", a matrix whose rows are the IDs, and
# "row_data", a data.frame with the descriptions.
.v2_read_functional_table <- function(path){
    first <- readLines(path, n = 1L)
    if( length(first) == 0L ){
        return(NULL)
    }
    header <- strsplit(first, "\t", fixed = TRUE)
    header <- header[[1L]]
    id_col <- max(match("description", header) - 1L, 1L, na.rm = TRUE)
    # The ID, description and category are text, and the others are counts
    is_annotation <- header %in% c("description", "category")
    is_annotation[[id_col]] <- TRUE
    col_classes <- ifelse(is_annotation, "character", "numeric")
    tab <- read.delim(
        path, check.names = FALSE, comment.char = "", quote = "",
        colClasses = col_classes)
    # An empty table has no data
    if( nrow(tab) == 0L ){
        return(NULL)
    }
    ids <- tab[[id_col]]
    mat <- tab[, !is_annotation, drop = FALSE] |> as.matrix()
    # Keep duplicated column names, like in .v2_read_table()
    colnames(mat) <- header[!is_annotation]
    mat[is.na(mat)] <- 0
    rownames(mat) <- ids
    row_data <- .v2_functional_row_data(
        ids, tab[["description"]], tab[["category"]])
    res <- list(counts = mat, row_data = row_data)
    return(res)
}

# Read the functional table of a single analysis. The tables have different
# formats:
# - pipeline v6: header, ID, description, GO category in GO tables, and
#   count, e.g.,
#     go          term                category            count
#     GO:0016020  membrane            cellular_component  2277
# - pipelines v3-v5: no header, and either count, ID and description, e.g.,
#     "55","IPR004995","Bacillus/Clostridium Ger spore germination protein"
#   or, in GO tables, ID, description, GO category and count, e.g.,
#     "GO:0009847","spore germination","biological_process","55"
# Returns a list like .v2_read_functional_table() with a one-column matrix.
.v2_read_functional_analysis_table <- function(path){
    # The tables can be gzipped. readLines() decompresses them.
    first <- readLines(path, n = 1L, warn = FALSE)
    # An empty table has no data
    if( length(first) == 0L ){
        return(NULL)
    }
    if( grepl("\t", first, fixed = TRUE) ){
        # Pipeline v6
        tab <- read.delim(
            path, check.names = FALSE, quote = "", comment.char = "",
            colClasses = "character")
        fields <- list(
            ids = tab[[1L]], description = tab[[2L]],
            category = tab[["category"]], counts = tab[["count"]])
    } else{
        fields <- .v2_read_quoted_functional_table(path)
    }
    if( length(fields[["counts"]]) == 0L ){
        return(NULL)
    }
    counts <- .v2_as_numeric(fields[["counts"]])
    if( all(is.na(counts)) ){
        stop("The table has no counts.", call. = FALSE)
    }
    counts[is.na(counts)] <- 0
    mat <- matrix(counts, dimnames = list(fields[["ids"]], NULL))
    row_data <- .v2_functional_row_data(
        fields[["ids"]], fields[["description"]], fields[["category"]])
    res <- list(counts = mat, row_data = row_data)
    return(res)
}

# Read the functional table of an analysis from pipelines v3-v5. The fields
# are quoted and separated by commas. Descriptions can have quotes that are
# not escaped, e.g., "streptomycin 3"-kinase", so read.csv() cannot be used.
# Instead, the lines are split by '","'. Returns a list with the fields "ids",
# "description", "category" (NULL if not a GO table) and "counts".
.v2_read_quoted_functional_table <- function(path){
    lines <- readLines(path, warn = FALSE)
    lines <- lines[nzchar(lines)]
    lines <- sub('^"', "", lines)
    lines <- sub('"$', "", lines)
    rows <- strsplit(lines, '","', fixed = TRUE)
    field <- function(i) vapply(rows, function(row) row[i], character(1))
    # GO tables have ID, description, GO category and count, and the other
    # tables count, ID and description
    if( length(lines) > 0L && startsWith(lines[[1L]], "GO:") ){
        res <- list(
            ids = field(1L), description = field(2L), category = field(3L),
            counts = field(4L))
    } else{
        res <- list(
            ids = field(2L), description = field(3L), category = NULL,
            counts = field(1L))
    }
    return(res)
}

# Descriptions of the rows of a functional table, e.g., the names of GO terms.
# GO categories are written as "biological process" or "biological_process",
# so they are written the same way. An ID can be on several rows with
# different descriptions, e.g., in InterPro tables of pipeline v1, where
# "IPR000008" is both "C2 domain" and "C2 calcium-dependent membrane
# targeting", as the entry was renamed between analyses. The counts of such
# rows are summed when the profiles are combined (see .v2_build_tse()), and
# the first description is used, as when the descriptions of many tables are
# combined (see .v2_combine_row_data()).
.v2_functional_row_data <- function(ids, description, category = NULL){
    category <- category %||% rep(NA_character_, length(ids))
    category <- gsub(" ", "_", category, fixed = TRUE)
    is_first <- !duplicated(ids)
    row_data <- data.frame(
        description = description[is_first], category = category[is_first],
        row.names = ids[is_first])
    return(row_data)
}

# Combine the descriptions of the rows of many tables. The same ID can be in
# many tables, and the description from the first table is used.
.v2_combine_row_data <- function(row_data){
    combined <- do.call(rbind, unname(row_data))
    ids <- lapply(row_data, rownames) |> unlist()
    combined <- combined[!duplicated(ids), , drop = FALSE]
    rownames(combined) <- ids[!duplicated(ids)]
    # Drop the category if the table has none, e.g., in InterPro tables
    if( is.na(combined[["category"]]) |> all() ){
        combined[["category"]] <- NULL
    }
    return(combined)
}

# Lineages are written differently in different tables, e.g.,
# "Root;k__Bacteria;p__Firmicutes" and "k__Bacteria; p__Firmicutes". Use the
# same format, so that the same taxa can be matched between tables.
.v2_clean_lineage <- function(lineage){
    lineage <- trimws(lineage)
    # Remove spaces after ";"
    lineage <- gsub(";\\s+", ";", lineage)
    lineage <- sub("^Root;", "", lineage)
    return(lineage)
}

# Taxonomy ranks by their prefixes in lineages, from the highest to the lowest
# rank. PR2 has its own ranks: supergroup, division and subdivision, e.g.,
# "d__Bacteria;sg__Terrabacteria;dv__Firmicutes;sdv__Firmicutes_X;c__Bacilli".
.V2_RANKS <- c(
    d = "domain", sk = "superkingdom", k = "kingdom", sg = "supergroup",
    dv = "division", sdv = "subdivision", p = "phylum", c = "class",
    o = "order", f = "family", g = "genus", s = "species")

# Parse taxonomy lineages into a table of taxonomy ranks. For example,
# "sk__Bacteria;k__;p__Bacillota" gives superkingdom "Bacteria" and phylum
# "Bacillota". Empty levels, e.g., "k__", and levels without a rank prefix,
# e.g., "Root", are skipped. ranks are the ranks that the result always has,
# e.g., the ranks of the reference database, so that data of different
# studies has the same columns. Other ranks are included if any taxon has
# them.
#' @importFrom S4Vectors DataFrame
.v2_parse_lineage <- function(lineage, ranks = NULL){
    all_ranks <- .V2_RANKS
    mat <- matrix(
        NA_character_, nrow = length(lineage), ncol = length(all_ranks),
        dimnames = list(lineage, all_ranks))
    for( prefix in names(all_ranks) ){
        # Find the level of the rank, e.g., "p__Bacillota" for prefix "p".
        # The level is at the start of the lineage or after ";", so that,
        # e.g., "k__" does not match "sk__".
        pattern <- paste0("(^|;)", prefix, "__([^;]+)")
        matches <- regexec(pattern, lineage)
        hits <- regmatches(lineage, matches)
        # Each hit has the whole match, the separator and the name
        taxa <- vapply(hits, function(hit){
            if( length(hit) == 0L ){
                return(NA_character_)
            }
            return(trimws(hit[[3L]]))
        }, character(1))
        mat[, all_ranks[[prefix]]] <- taxa
    }
    # Keep the given ranks and the ranks that any taxon has
    is_found <- !is.na(mat)
    has_rank <- colSums(is_found) > 0L | all_ranks %in% ranks
    mat <- mat[, has_rank, drop = FALSE]
    res <- DataFrame(mat, check.names = FALSE)
    return(res)
}

# Combine abundance tables into a single TreeSE. tables are matrices whose
# rows are lineages and columns are analyses. col_data has the analyses as
# rows. Taxa are matched between the tables by their lineage, and a taxon
# that is missing from a table has zero count in its analyses. row_data
# describes the rows of functional tables, whose rows are IDs, e.g., GO
# terms. If it is NULL, rowData has the taxonomy ranks from the lineages, and
# ranks are the ranks that rowData always has (see .v2_parse_lineage()).
#' @importFrom Matrix sparseMatrix
#' @importFrom S4Vectors SimpleList DataFrame
.v2_build_tse <- function(
        tables, col_data, max.dense = 1e7, row_data = NULL, ranks = NULL){
    # Rows of the result are all the taxa of all the tables, and columns are
    # the analyses
    lineages <- lapply(tables, rownames) |> unlist() |> as.character() |>
        unique()
    features <- .v2_unique_taxa(lineages)
    feature_keys <- .v2_taxon_key(features)
    analyses <- rownames(col_data)

    # Most counts are zero, so the matrix is built as sparse to save memory.
    # A sparse matrix is built from its non-zero values, each given as
    # (row, column, count). For example, (2, 3, 5) means that taxon 2 has
    # count 5 in analysis 3.
    nonzero <- lapply(tables, function(tab){
        # Positions of non-zero values in this table as (row, column) pairs
        idx <- which(tab != 0, arr.ind = TRUE)
        taxa <- rownames(tab)[idx[, 1L]] |> .v2_taxon_key()
        ids <- colnames(tab)[idx[, 2L]]
        # Convert them to positions in the result. as.integer() and
        # as.numeric() are needed for tables without non-zero values.
        res <- data.frame(
            row = match(taxa, feature_keys) |> as.integer(),
            col = match(ids, analyses) |> as.integer(),
            count = as.numeric(tab[idx]))
        return(res)
    })
    nonzero <- do.call(rbind, nonzero)
    # sparseMatrix() sums the counts of the same (row, column) pair. This
    # sums the rows of a table that are the same taxon or ID, e.g., an ID on
    # several rows of a functional table (see .v2_functional_row_data()) or
    # lineages that are the same after cleaning (see .v2_clean_lineage()).
    counts <- sparseMatrix(
        i = nonzero[["row"]], j = nonzero[["col"]], x = nonzero[["count"]],
        dims = c(length(features), length(analyses)),
        dimnames = list(features, analyses))
    # Dense matrices are more widely supported, so use one unless the data
    # is large
    if( prod(dim(counts)) <= max.dense ){
        counts <- as.matrix(counts)
    }

    assays <- SimpleList(counts = counts)
    if( is.null(row_data) ){
        row_data <- .v2_parse_lineage(features, ranks)
    } else{
        row_data <- row_data[features, , drop = FALSE] |>
            DataFrame(check.names = FALSE)
    }
    col_data <- .v2_drop_empty_metadata(col_data) |>
        DataFrame(check.names = FALSE)
    tse <- TreeSummarizedExperiment(
        assays = assays, rowData = row_data, colData = col_data)
    return(tse)
}

# The same taxon can be written with or without spaces in different tables,
# e.g., "g__Candidatus Nitrososphaera" in the tables of analyses and
# "g__CandidatusNitrososphaera" in the summary tables of pipelines v2-3. Taxa
# are therefore matched by their lineage without spaces.
.v2_taxon_key <- function(lineage){
    key <- gsub(" ", "", lineage, fixed = TRUE)
    return(key)
}

# Get the unique taxa of lineages, matched by .v2_taxon_key(). Each taxon is
# named by the spelling that has the most spaces, e.g.,
# "g__Candidatus Nitrososphaera", because it is the easiest to read. The taxa
# are in the order in which they first appear.
.v2_unique_taxa <- function(lineages){
    keys <- .v2_taxon_key(lineages)
    # Sort the spellings so that the one with the most spaces is the first
    n_spaces <- nchar(lineages) - nchar(keys)
    best_first <- order(n_spaces, decreasing = TRUE)
    is_best <- !duplicated(keys[best_first])
    best <- lineages[best_first][is_best]
    names(best) <- keys[best_first][is_best]
    # Keep the order of first appearance
    taxa <- best[unique(keys)] |> unname()
    return(taxa)
}
