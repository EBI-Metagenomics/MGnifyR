# Internal function to actually perform the http request. Build up the URL then issues
# a GET, parsing the returned JSON into a nested list (uses \code{jsonlite} internally?)
# Previously cached results may be retrieved from disk without resorting to calling the MGnify server.

#' Low level MGnify API handler
#'
#' \code{mgnify_retrieve_json} deals with handles the actual HTTP GET calls for the MGnifyR package, handling API pagination,
#' local result caching, and authentication cookies for access
#' to restricted or pre-release datasets.Although principally intended for internal MGnifyR use , it's exported for direct invocation.
#' Generally though it's not recommended for use by users.
#'
#' @importFrom urltools parameters
#' @importFrom httr add_headers
#' @importFrom httr GET
#' @importFrom httr config
#' @importFrom httr content
#'
#' @param client MGnifyR client
#' @param path top level search point for the query. One of \code{biomes}, \code{samples}, \code{runs} etc. Basically includes
#' all parts of the URL between the base API url and the parameter specifications
#' @param complete_url \emph{complete} url to search, usuaally retrieved from previous query in the "related" section.
#' @param qopts named list or vector containing options/filters to be URL encoded and appended to query as key/value pairs
#' @param maxhits Maxmium number of data entries to return. The actual number of hits returned may be higher than this value,
#' as this parameter only clamps after each full page is processed. Set to <=0 to disable - i.e. retrieve all items.
#' @param usecache Should successful queries be cached on disk locally? There are unresolved questions about whether this is
#' a sensible thing to do, but it remains as an option. It probably makes sense for single accession grabs, but not for
#' (filtered) queries - which are liable to change as new data is added to MGnify. Also caching only works for the first page.
#' @param Debug Should we print out lots of information while doing the grabbing?
#' @return \code{list} of results after pagination is dealt with.
#' @export
mgnify_retrieve_json <- function(client, path="biomes", complete_url=NULL, qopts=NULL,maxhits=200, usecache = F, Debug=F){


    #client@warnings turns on debugging too:

    if(client@warnings){
        Debug <- T
    }
    # Set up the base url
    # Are we using internal paths?
    if (is.null(complete_url)){
        fullurl <- paste(client@url, path, sep="/")
    }
    #Or direct links from e.g. a "related" section
    else{
        #Set the full url, but clean off any existing parameters (page, format etc) as they'll be added back later:
        fullurl <- complete_url
        parameters(fullurl) <- NULL
        path <- substr(fullurl, nchar(client@url) + 2, nchar(fullurl))
    }

    #cat(fullurl)

    #convert to csv if filters are lists.
    #This doesn't check if they ~can~ be searched for in the API,
    #which is an issue since no error is returned by the JSON if the search
    #is invalid - we only get a result as if no query was present...
    tmpqopts <- lapply(qopts,function(x) paste(x,collapse = ','))

    #Include the json and page position options
    #full_qopts <- as.list(c(format="json", tmpqopts, page=1))
    full_qopts <- as.list(c(format="json", tmpqopts))
    #Build up the cache name anyway - even if it's not ultimately used:
    fname_list <- c(path, names(unlist(full_qopts)), unlist(full_qopts))
    cache_fname <- paste(fname_list,collapse = "_")
    cache_full_fname <- paste(client@cache_dir, '/', cache_fname, '.RDS', sep="")


    ## Quick check to see if we should clear the disk cache ~for this specific call~ - used for debugging
    # and when MGnify breaks
    if(usecache & client@clear_cache){
        message(paste("clear_cache is TRUE: deleting ", cache_full_fname, sep=""))
        tryCatch(unlink(cache_full_fname), error=warning)
    }

    # Do we want to try and use a cache to speed things up?
    if(usecache & file.exists(cache_full_fname)){
            final_data <- readRDS(cache_full_fname)
    }else{

        #Authorization: Bearer <your_token>
        if(!is.null(client@authtok)){
            add_headers(.headers = c(Authorization = paste("Bearer", client@authtok, sep=" ")))
        }
        res <- GET(url=fullurl, config(verbose=Debug), query=full_qopts )
        data <- content(res)

        #At this point, data$data is either a list of lists or a single named list. If it's a single entry, it needs embedding in
        #a list for consistency downstream
        #datlist is built up as a list of pages, where each entry must be another list. Thus, on the first page,
        #
        datlist <- list()
        if (!is.null(names(data$data))){
        #Create something to store the returned data

            datlist[[1]] <- list(data$data)
        }else{
            datlist[[1]] <- data$data
        }
            #cat(str(data))
        # Check to see if there's pagination required
        if ("meta" %in% names(data)){
            #Yes, paginate
            pstart <- as.numeric(data$meta$pagination$page)
            pend <- as.numeric(data$meta$pagination$pages)

            for (p in seq(pstart+1,pend)){    # We've already got the first one

                full_qopts$page <- p
                if(!is.null(client@authtok)){
                    add_headers(.headers = c(Authorization = paste("Bearer", client@authtok, sep=" ")))
                }
                curd <- content(GET(fullurl, config(verbose=Debug), query=full_qopts ))
                datlist[[p]] <- curd$data
                #Check to see if we've pulled enough entries
                if(maxhits > 0){
                    curlen <- sum(sapply(datlist, length))
                    if (curlen > maxhits){
                        break
                    }
                }
            }
        }
        #if(length(datlist) > 1){
        final_data <- unlist(datlist, recursive=F)

        if (usecache && !file.exists(cache_full_fname)){
            #Make sure the directory is created...
        dir.create(dirname(cache_full_fname), recursive = T, showWarnings = client@warnings)
            saveRDS(final_data, file = cache_full_fname)
        }
    }
    final_data
}
