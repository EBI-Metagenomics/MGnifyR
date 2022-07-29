# Constructor to allow logging in with username/password
#' Instantiate the MGnifyR client object
#'
#' All functions in the MGnifyR package take a \code{mgnify_client} object as their first argument. While not essential
#' to querying the raw MGnify API (which is exposed as relative standard JSONAPI), the object allows the simple handling of both
#' user authentication and access to private data, and local on-disk caching of results.
#'
#' @importFrom httr POST
#' @importFrom httr content
#'
#' @param url (To be described)
#' @param username optional username to authenticate.
#' @param password optional password for authentication.
#' @param usecache whether to enable on-disk caching of results during this session. In most use cases should be TRUE.
#' @param cache_dir specifies a folder to contain the local cache. If NULL, and usecache is TRUE, the new subdirectory \code{.MGnifyR_cache}
#' in the current working directory will be used. Note that cached files are persistent, so the cache directory may be reused between sessions,
#' taking advantage of previously downloaded results. The directory will be created if it doesn't exist already.
#' @param warnings debug flag to print extra output during invocation of some MGnifyR functions. Defaults to FALSE.
#' @param use_memcache flag to indicate whether functional results obtained when \code{bulk_dl} is \code{TRUE} in \code{mgnify_get_analyses_results} should
#' be stored in an in-memory cache, rather than the cached input being re-read for each accession. this is currently NOT working
#' properly and should therefore be set \code{FALSE} (the default). It has the potential to speed up searches considerably though, especially
#' for studies with a large number of samples, so will be implemented properly in the future.
#' @examples
#' my_client <- mgnify_client(username="Webin-1122334", password="SecretPassword", usecache=T, cache_dir = "/scratch/MGnify_cache_location")
#' @export
mgnify_client <- function(url=NULL,username=NULL,password=NULL,usecache=F,cache_dir=NULL, warnings=F, use_memcache=F){
    if (is.null(url)){
        url <- baseurl
    }

    authtok <- NA_character_

    #Check to see if we're goint to try and get an authentication token:
    if (!is.null(username) && !is.null(password)){
        r <- httr::POST(paste(url, "utils/token/obtain", sep="/"),
                                     body=list(username=username, password=password),
                                     encode="json")
        cont <- httr::content(r)
        if ("data" %in% names(cont)){
            authtok <- cont$data$token
        }
        else{
            stop("Failed to authenticate")
        }
    }
    #Assume we're not using it
    cachepath <- NA_character_
    if(usecache){
        if (is.null(cache_dir) ){
            cachepath <- paste(getwd(),'.MGnifyR_cache',sep="/")
        }else{
            cachepath <- cache_dir
        }
        #Make it if needed - assume the user is sensible and the path will work...
        dir.create(cachepath,showWarnings = F)
    }

    #Return the final object
    #@importFrom methods new
    new("mgnify_client", url=url, authtok=authtok, cache_dir = cachepath, warnings=warnings, memcache=list(), use_memcache=use_memcache)
}
