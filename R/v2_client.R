# Access to MGnify API v2: authentication and downloading with caching

# Base URL of the MGnify API v2
.V2_URL <- "https://www.ebi.ac.uk/metagenomics/api/v2/"
# BioSamples has all the attributes of samples, e.g., "disease status"
.BIOSAMPLES_BULK_URL <- paste0(
    "https://www.ebi.ac.uk/biosamples/v2/samples/bulk-fetch?accessions=")

# Session state, i.e., authentication token for private data
.v2_state <- new.env(parent = emptyenv())

#' Log in to MGnify to access private data
#'
#' @details
#' Most MGnify data is public and can be fetched without logging in. To access
#' private data, log in with your ENA Webin credentials. The login lasts for
#' the R session, and it is used by all the functions that fetch data from
#' MGnify API v2, such as \code{\link{searchStudies}} and
#' \code{\link{fetchStudies}}.
#'
#' @param username A single character value specifying the Webin username,
#' e.g., \code{"Webin-1122334"}.
#'
#' @param password A single character value specifying the password.
#'
#' @return \code{TRUE} invisibly, if login was successful.
#'
#' @examples
#' # Wrong credentials give an error
#' try(mgnifyLogin("Webin-0", "wrong"))
#'
#' # With your own ENA Webin credentials, log in and fetch private data
#' # mgnifyLogin("Webin-1122334", "SecretPassword")
#' # tse <- fetchStudies("MGYS00000001")
#'
#' @export
#' @importFrom curl new_handle handle_setopt handle_setheaders
#' @importFrom curl curl_fetch_memory
#' @importFrom jsonlite toJSON parse_json
mgnifyLogin <- function(username, password){
    .check_string(username, "username")
    .check_string(password, "password")
    handle <- new_handle()
    handle_setheaders(handle, "Content-Type" = "application/json")
    # Send the credentials as JSON, e.g., {"username": ..., "password": ...}
    credentials <- list(username = username, password = password) |>
        toJSON(auto_unbox = TRUE)
    handle_setopt(handle, postfields = credentials)
    url <- paste0(.V2_URL, "auth/sliding")
    res <- curl_fetch_memory(url, handle)
    # The response has the token if the login was successful
    token <- tryCatch({
        response <- rawToChar(res[["content"]]) |> parse_json()
        response[["token"]]
    }, error = function(e) NULL)
    if( is.null(token) ){
        stop("Failed to log in.", call. = FALSE)
    }
    .v2_state[["token"]] <- token
    return(invisible(TRUE))
}

# Downloaded files are cached for the R session. To keep the cache between
# sessions, the user can set options(MGnifyR.cache = "path/to/cache").
.v2_cache_dir <- function(){
    default_dir <- file.path(tempdir(), "MGnifyR")
    cache_dir <- getOption("MGnifyR.cache", default_dir)
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    return(cache_dir)
}

# Download URLs in parallel through the cache. Returns local file paths in the
# same order as urls, NA for those that could not be downloaded.
#' @importFrom curl multi_download
#' @importFrom tools md5sum
.v2_download <- function(urls, ext){
    if( length(urls) == 0L ){
        return(character(0))
    }
    # Each URL has its own file in the cache, named by the MD5 hash of the
    # URL, e.g., "4e8a5ff852cbd0798d914ae7bc2d1691.tsv". The URL itself cannot
    # be used as the file name: it has characters that are not allowed in file
    # names, e.g., "/", ":" and "?", and it can be too long for Windows. The
    # last part of the URL is not unique, e.g., the same file name is used in
    # many studies. The hash is unique, short and always a valid file name.
    hashes <- vapply(urls, function(url){
        bytes <- charToRaw(url)
        hash <- md5sum(bytes = bytes)
        return(hash)
    }, character(1))
    files <- file.path(.v2_cache_dir(), paste0(hashes, ext))
    # Download only the files that are not in the cache yet
    in_cache <- file.exists(files)
    todo <- urls[!in_cache] |> unique()
    dest <- files[match(todo, urls)]
    .v2_download_retry(todo, dest)
    files[!file.exists(files)] <- NA_character_
    return(files)
}

# Download URLs in parallel to files. Try up to 5 times, because downloads can
# fail temporarily, e.g., when the FTP server refuses connections for a
# while. Waiting is cheap compared to a missing file: if a summary table of a
# study cannot be downloaded, its analyses are fetched one by one, which can
# take minutes. Files that could not be downloaded are not created.
.v2_download_retry <- function(urls, files){
    for( attempt in seq_len(5L) ){
        if( length(urls) == 0L ){
            break
        }
        # Wait before trying again, longer each time: 2, 4, 8 and 16 seconds,
        # i.e., at most 30 seconds in total
        if( attempt > 1L ){
            Sys.sleep(2^(attempt - 1L))
        }
        # Download to temporary ".part" files, and rename them only after a
        # successful download, so that broken downloads are not cached
        tmp <- paste0(files, ".part")
        res <- .v2_multi_download(urls, tmp)
        # "success" means only that the transfer completed, so check also the
        # HTTP status, e.g., 200 is OK and 404 is "not found"
        ok <- res[["success"]] %in% TRUE & res[["status_code"]] %in% 200:299
        file.rename(tmp[ok], files[ok])
        unlink(tmp[!ok])
        # Do not try again if the file does not exist or is not accessible
        # (HTTP status 4xx). Status 429 "too many requests" is tried again.
        permanent <- res[["status_code"]] %in% c(400:428, 430:499)
        retry <- !ok & !permanent
        urls <- urls[retry]
        files <- files[retry]
    }
    return(invisible(NULL))
}

# Timeouts of downloads in seconds. Without them, a connection that the
# server never answers, or a transfer that stops, would wait forever and the
# whole fetch would hang. A connection must be made in 30 seconds, and a
# transfer that gets no data for 60 seconds is stopped and tried again (see
# .v2_download_retry()). Large files are not limited, as long as data comes.
.V2_CONNECT_TIMEOUT <- 30L
.V2_STALL_TIMEOUT <- 60L

# Download URLs in parallel to files. If the user has logged in, the login
# token is sent to MGnify API to get private data, but not to other servers,
# e.g., BioSamples or the FTP server of MGnify. Returns a data.frame with the
# columns "success" and "status_code" of each URL.
.v2_multi_download <- function(urls, files){
    # Show a progress bar for large downloads
    progress <- interactive() && length(urls) > 10L
    # multi_download() sends the same header to all URLs, so the URLs of
    # MGnify API and other servers are downloaded separately
    is_api <- startsWith(urls, .V2_URL)
    res <- data.frame(
        success = rep(NA, length(urls)),
        status_code = rep(NA_real_, length(urls)))
    # Other servers, without the login token
    if( any(!is_api) ){
        downloaded <- multi_download(
            urls[!is_api], files[!is_api], progress = progress,
            connecttimeout = .V2_CONNECT_TIMEOUT,
            low_speed_limit = 1L, low_speed_time = .V2_STALL_TIMEOUT)
        res[!is_api, ] <- downloaded[, c("success", "status_code")]
    }
    # MGnify API, with the login token if the user has logged in. Without
    # the token, the header is empty.
    token <- .v2_state[["token"]]
    header <- character(0)
    if( !is.null(token) ){
        header <- paste("Authorization: Bearer", token)
    }
    if( any(is_api) ){
        downloaded <- multi_download(
            urls[is_api], files[is_api], progress = progress,
            httpheader = header, connecttimeout = .V2_CONNECT_TIMEOUT,
            low_speed_limit = 1L, low_speed_time = .V2_STALL_TIMEOUT)
        res[is_api, ] <- downloaded[, c("success", "status_code")]
    }
    return(res)
}
