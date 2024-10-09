result <- searchAnalysis(mg, "studies", c("MGYS00005058"))

metadata <- getMetadata(mg, result)
samples <- metadata$sample_biosample[1:2]

# WHat files are available?
?accession=SAMEA10104908&result=read_run&format=json

# Download files.

searchEnaFile <- function(accession, wide.format = FALSE, ...){
    res <- lapply(accession, function(x){
        .search_ena_files(accession = x, ...)
    })
    res <- do.call(rbind, res) |> as.data.frame()
    
    if( wide.format ){
        res <- res |>
            group_by(run_accession, accession, file_type) |>
            mutate(row_id = row_number()) |>
            ungroup() |>
            pivot_wider(
                names_from = c(file_type, row_id),  # Create unique columns for each file_type based on row_id
                values_from = file_name,
                names_sep = "_"  # Add a separator for clarity
            ) |> as.data.frame()
    }
    return(res)
}

# values as list so that user can have more options than these?? type is only needed (and format).
.search_ena_files <- function(
        accession,
        ena.url = "https://www.ebi.ac.uk/ena/portal/api/",
        type = "filereport",
        result = "read_run",
        format = "json",
        ...){
    url <- paste0(
        ena.url,
        type, "?",
        "result=", result, "&",
        "format=", format, "&",
        "accession=", accession
        )
    res <- GET(url = url)
    # Get the data
    res <- content(res) |> unlist()
    if( !is.null(res) ){
        res <- lapply(names(res), function(x){
            temp <- strsplit(res[[x]], ";") |> unlist()
            names(temp) <- rep(x, length(temp))
            return(temp)
        })
        res <- res |> unlist()
        res[["accession"]] <- accession
        colnames <- names(res)
        res <- as.data.frame(res) |> t() |> as.data.frame()
        colnames(res) <- colnames
        res <- res |> pivot_longer(
            cols = -c(accession, run_accession), # Exclude accession columns from pivoting
            names_to = "file_type",
            values_to = "file_name"
        )
    }
    
    return(res)
}

files <- searchEnaFile(c("SAMEA10104908", "SAMEA10104909"))
files <- searchEnaFile(c("SAMEA10104908", "SAMEA10104909"), wide.format = TRUE)
url <- files[[1]]
getFile(mg, url)


