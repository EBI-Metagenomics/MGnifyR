#' Constructor for creating a MgnifyClient object to allow the access to
#' MGnify database.
#'
#' @details
#' All functions in the MGnifyR package take a \code{MgnifyClient} object as
#' their first argument. While not essential to querying the raw MGnify API
#' (which is exposed as relative standard JSONAPI), the object allows the
#' simple handling of both user authentication and access to private data,
#' and local on-disk caching of results.
#'
#' @param username A single character value specifying an optional username for
#' authentication. (By default: \code{username = NULL})
#'
#' @param password A single character value specifying an optional password for
#' authentication. (By default: \code{password = NULL})
#'
#' @param useCache A single boolean value specifying whether to enable on-disk
#' caching of results during this session. In most use cases should be TRUE.
#' (By default: \code{useCache = FALSE})
#'
#' @param cacheDir A single character value specifying a folder to contain the
#' local cache. Note that cached files are persistent, so the cache directory
#' may be reused between sessions, taking advantage of previously downloaded
#' results. The directory will be created if it doesn't exist already.
#' (By default: \code{cacheDir = tempdir()})
#'
#' @param showWarnings A single boolean value specifying whether to print
#' warnings during invocation of some MGnifyR functions.
#' (By default: \code{showWarnings = FALSE})
#'
#' @param verbose A single boolean value specifying whether to print extra
#' output during invocation of some MGnifyR functions.
#' (By default: \code{verbose = FALSE})
#'
#' @param clearCache A single boolean value specifying whether to clear the
#' cache. (By default: \code{clearCache = FALSE})
#'
#' @param ... optional arguments:
#' \itemize{
#'   \item \strong{url} A single character value specifying an url address of
#'   the database. (By default:
#'   \code{url = "https://www.ebi.ac.uk/metagenomics/api/v1"})
#' }
#'
#' @return A MgnifyClient object.
#'
#' @examples
#' my_client <- MgnifyClient(
#'     useCache = TRUE, cacheDir = "/scratch/MGnify_cache_location"
#'     )
#'
#' \dontrun{
#' # Use username and password to get access to non-public data
#' my_client <- MgnifyClient(
#'     username = "Webin-1122334", password = "SecretPassword",
#'     useCache = TRUE, cacheDir = "/scratch/MGnify_cache_location"
#'     )
#'}
#'
#' @name MgnifyClient
NULL

#' @rdname MgnifyClient
#' @importFrom methods new
#' @export
MgnifyClient <- function(
        username = NULL, password = NULL, useCache = FALSE,
        cacheDir = tempdir(), showWarnings = FALSE, verbose = TRUE,
        clearCache = FALSE, ...){
    ############################### INPUT CHECK ################################
    if( !(is.null(username) || .is_non_empty_string(username)) ){
        stop(
            "'username' must be NULL or single character value specifying ",
            "the username.", call. = FALSE)
    }
    if( !(is.null(password) || .is_non_empty_string(password)) ){
        stop(
            "'password' must be NULL or single character value specifying ",
            "the password.", call. = FALSE)
    }
    if( !.is_a_bool(useCache) ){
        stop(
            "'useCache' must be a boolean value specifying whether to use ",
            "on-disk caching.", call. = FALSE)
    }
    if( !.is_non_empty_string(cacheDir) ){
        stop(
            "'cacheDir' must be single character value specifying ",
            "the the directory for cache.", call. = FALSE)
    }
    if( !.is_a_bool(showWarnings) ){
        stop(
            "'showWarnings' must be a boolean value specifying whether print ",
            "warnings during invocation of MGnifyR functions.",
            call. = FALSE)
    }
    if( !.is_a_bool(verbose) ){
        stop(
            "'verbose' must be a boolean value specifying whether print ",
            "extra output during invocation of MGnifyR functions.",
            call. = FALSE)
    }
    if( !.is_a_bool(clearCache) ){
        stop(
            "'clearCache' must be a boolean value specifying whether to ",
            "clear the cache.", call. = FALSE)
    }
    ############################# INPUT CHECK END ##############################
    # Get the url address
    url <- .get_url_address(...)
    # Authentication token is NA as default
    authtok <- NA_character_
    # Check to see if we're going to try and get an authentication token:
    if (!is.null(username) && !is.null(password)){
        # Fetch username vs password data from database
        r <- POST(
            paste(url, "utils/token/obtain", sep = "/"),
            body = list(username = username, password = password),
            encode = "json")
        # If the authentication was not successful, returned value do not
        # include data
        cont <- content(r, ...)
        if ("data" %in% names(cont)){
            authtok <- cont$data$token
        } else{
            stop("Failed to authenticate.", call. = FALSE)
        }
    }
    # Get the directory where cache will be stored.
    # If user has specified the subdirectory, ensure that it works in any
    # system by adding correct "/".
    cacheDir <- as.list(strsplit(cacheDir, "[/\\\\]")[[1]])
    cacheDir <- do.call(file.path, cacheDir)
    # Add subdirectory. If user has specified for example working directory,
    # the directory would be full of files. This is unintentional.
    cacheDir <- file.path(cacheDir, ".MGnifyR_cache")
    # Make it if needed - assume the user is sensible and the path will
    # work...
    if( useCache ){
        dir.create(cacheDir, recursive = TRUE, showWarnings = FALSE)
    }
    # Return the final object
    obj <- new(
        "MgnifyClient",
        databaseUrl = url,
        authTok = authtok,
        useCache = useCache,
        cacheDir = cacheDir,
        showWarnings = showWarnings,
        clearCache = clearCache,
        verbose = verbose
    )
    return(obj)
}

################################ HELP FUNCTIONS ################################

# This function is just to remove url from main function's arguments.
.get_url_address <- function(
        url = "https://www.ebi.ac.uk/metagenomics/api/v1", ...){
    ############################### INPUT CHECK ################################
    if( !(.is_non_empty_string(url)) ){
        stop(
            "'url' must be a single character value specifying ",
            "the URL address.", call. = FALSE)
    }
    ############################# INPUT CHECK END ##############################
    return(url)
}

