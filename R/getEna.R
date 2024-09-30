result <- searchAnalysis(mg, "studies", c("MGYS00005058"))

metadata <- getMetadata(mg, result)
samples <- metadata$sample_biosample[1:2]

# WHat files are available?
?accession=SAMEA10104908&result=read_run&format=json

# Download files.

searchEnaFile <- function(accession, ...){
    res <- lapply(accession, function(x){
        .search_ena_files(accession = x, ...)
    })
    accession <- accession[ lengths(res) > 0 ]
    res <- dplyr::bind_rows(res)
    res[["accession"]] <- accession
    res <- apply(res, 1, function(row) unlist(strsplit(row, ";")))
    res <- res |> t() |> as.data.frame()
    return(res)
}

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
    file_data <- content(res) |> unlist()
    # file_data <- lapply(file_data, function(x){
    #     strsplit(x, ";")
    # })
    # file_data <- file_data |> unlist()
    return(file_data)
}

files <- searchEnaFile(c("SAMEA10104908", "SAMEA10104909"))
url <- files[[1]]
getFile(mg, url)


