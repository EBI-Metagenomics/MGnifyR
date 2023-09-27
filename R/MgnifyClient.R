#' A MgnifyClient object
#'
#' @details An object that are required by functions of MGnifyR package.
#'
#' @slot url A single character value specifying an URL address of database.
#'
#' @slot authTok A single character value specifying authentication token.
#'
#' @slot cacheDir A single character value specifying cache directory.
#'
#' @slot warnings A single boolean value specifying whether to show warnings.
#'
#' @slot useMemCache A single boolean value specifying whether to use on-disk
#' memory.
#'
#' @slot memCache A single character value specifying on-disk memory directory.
#'
#' @slot clearCache A single boolean value specifying whether to clear cache.
#'
#' @section Constructor: see \code{MgnifyClient}.
#'
#' @name MgnifyClient
NULL

#' @rdname MgnifyClient
#' @importFrom httr POST
#' @importFrom httr content
#' @export
MgnifyClient <- setClass(
    "MgnifyClient", slots = list(url = "character", authTok = "character",
                                  cacheDir = "character", warnings = "logical",
                                  useMemCache = "logical", memCache = "list",
                                  clearCache = "logical"),
    prototype = list(url = "https://www.ebi.ac.uk/metagenomics/api/v1",
                     authTok = NULL, cacheDir = NULL,
                     useMemCache = FALSE, memCache = list(),
                     clearCache = FALSE))

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
#' local cache. If NULL, and useCache is TRUE, the new subdirectory
#' \code{.MGnifyR_cache} in the current working directory will be used. Note
#' that cached files are persistent, so the cache directory may be reused
#' between sessions, taking advantage of previously downloaded results. The
#' directory will be created if it doesn't exist already.
#' (By default: \code{cachedir = NULL})
#'
#' @param warnings A single boolean value specifying whether to print extra
#' output during invocation of some MGnifyR functions.
#' (By default: \code{warnings = FALSE})
#'
#' @param useMemCache A single boolean value specifying whether to indicate
#' whether functional results obtained when \code{bulk_dl} is \code{TRUE}
#' in \code{mgnify_get_analyses_results} should be stored in an in-memory
#' cache, rather than the cached input being re-read for each accession. this
#' is currently NOT working properly and should therefore be set \code{FALSE}.
#' It has the potential to speed up searches considerably though, especially
#' for studies with a large number of samples, so will be implemented properly
#' in the future. (By default: \code{useMemCache = FALSE})
#'
#' @param ... optional arguments:
#' \itemize{
#'   \item{url}{ A single character value specifying an url address of the
#'   the database.
#'   (By default: \code{url = "https://www.ebi.ac.uk/metagenomics/api/v1"})}
#' }
#'
#' @return A MgnifyClient object.
#'
#' @examples
#' \dontrun{
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
        username = NULL, password = NULL, useCache = FALSE, cacheDir = NULL,
        warnings = FALSE, useMemCache = FALSE, ...){
    ############################### INPUT CHECK ################################
    if( !(is.null(username) || .is_non_empty_string(username)) ){
        stop("'username' must be NULL or single character value specifying ",
             "the username.", call. = FALSE)
    }
    if( !(is.null(password) || .is_non_empty_string(password)) ){
        stop("'password' must be NULL or single character value specifying ",
             "the password.", call. = FALSE)
    }
    if( !.is_a_bool(useCache) ){
        stop("'useCache' must be a boolean value specifying whether to use ",
             "on-disk caching.", call. = FALSE)
    }
    if( !(is.null(cacheDir) || .is_non_empty_string(cacheDir)) ){
        stop("'cacheDir' must be NULL or single character value specifying ",
             "the the directory for cache.", call. = FALSE)
    }
    if( !.is_a_bool(warnings) ){
        stop("'wanings' must be a boolean value specifying whether print ",
             "extra output during invocation of MGnifyR functions.", call. = FALSE)
    }
    if( !.is_a_bool(useMemCache) ){
        stop("'useMemCache' must be a boolean value specifying whether use ",
             "on-disk memory.", call. = FALSE)
    }
    ############################# INPUT CHECK END ##############################
    # Get the url address
    url <- .get_url_address(...)
    # Authentication token is NA as default
    authtok <- NA_character_
    # Check to see if we're going to try and get an authentication token:
    if (!is.null(username) && !is.null(password)){
        # Fetch username vs password data from database
        r <- POST(paste(url, "utils/token/obtain", sep = "/"),
                  body = list(username=username,
                              password=password),
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
    # Assume we're not using it
    cachepath <- NA_character_
    # If user has specified that on-disk cache will be used
    if(useCache){
        if (is.null(cacheDir) ){
            cachepath <- paste(getwd(), ".MGnifyR_cache", sep = "/")
        } else{
            cachepath <- cacheDir
        }
        # Make it if needed - assume the user is sensible and the path will work...
        dir.create(cachepath, showWarnings = FALSE)
    }
    # Return the final object
    obj <- new("MgnifyClient", url = url, authTok = authtok,
               cacheDir = cachepath, warnings = warnings, memCache = list(),
               useMemCache = useMemCache)
    return(obj)
}

################################ HELP FUNCTIONS ################################

# This function is just to remove url from main function's arguments.
.get_url_address <- function(
        url = "https://www.ebi.ac.uk/metagenomics/api/v1", ...){
    ############################### INPUT CHECK ################################
    if( !(.is_non_empty_string(url)) ){
        stop("'url' must be a single character value specifying ",
             "the URL address.", call. = FALSE)
    }
    ############################# INPUT CHECK END ##############################
    return(url)
}

