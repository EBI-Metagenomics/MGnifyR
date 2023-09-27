################################### TESTING ###################################
# Methods for testing

.is_a_bool <- function(x){
    is.logical(x) && length(x) == 1L && !is.na(x)
}

.is_non_empty_character <- function(x){
    is.character(x) && all(nzchar(x))
}

.is_non_empty_string <- function(x){
    is.character(x) && length(x) == 1L
}

.is_an_integer <- function(x){
    is.numeric(x) && length(x) == 1L && x%%1==0
}

################################ HELP FUNCTIONS ################################
# Help functions that are utilized by multiple methods

########################## .mgnify_attr_list_to_df_row #########################
# Not exporting this - if people want to they can use the
# rjsonapi functionality. Internally, it takes the "attributes" list
# and converts it into a single row data.frame. For some entries, there is a
# sublist of key/value pairs. metadata_key allows these to be included as
# columns in the result.
.mgnify_attr_list_to_df_row <- function (json, metadata_key = NULL){
    # Get what kind of metadata the data includes
    attrlist <- names(json$attributes)
    # If the type of metadata is specified
    if (!is.null(metadata_key)){
        # Get metadata related to specific key
        metaattrlist <- json$attributes[[metadata_key]]
        metlist <- sapply(metaattrlist, function(x) x$value)
        names(metlist) <- sapply(metaattrlist, function(x) x$key)
        # Get metadata without the key
        baseattrlist <- attrlist[!(attrlist %in% c(metadata_key))]
        # Combine metadata
        df <- as.data.frame(t(unlist(c(
            json$attributes[baseattrlist], metlist))),
            stringsAsFactors = FALSE)
    }else{
        # Get all the metadata without key extraction
        df <- as.data.frame(t(unlist(json["attributes"])),
                            stringsAsFactors = FALSE)
    }
    # Add accession code and type of data
    df$accession <- json$id
    df$acc_type <- json$type
    # Add accession code also to rownames
    rownames(df) <- df$accession
    return(df)
}

############################## .mgnify_get_x_for_y #############################
# Helper function for getting relative paths in the API
# Not everything is implemented here - just what we
# need to get to the download or run areas
# Given an accession x, we want to get the link to get the url for the
# corresponding typeY JSONAPI path for child elements
#
# .mgnify_get_x_for_y determines the location of typeY child objects of x (typeX)
#
# This helper function, principally intended to be used internally,
# is used to match up related objects within the path. The inherently
# unhierarchical nature of the MGnify API makes it a bit inconsistent. This
# function acts as a quick way to determine how to get from one type to another,
# without having to special case within the code.
#
# Parameters:
# client MGnifyR client API object
# x Accession ID \code{char} of parent object
# typeX Type of accession \code{x}
# typeY Type of child object to return
# use.cache Whether to use on-disk cache
#
# Return:
# char complete url to access the result. Note this query is not run from here -
# just the URL is returned
#
# Examples:
# cl <- new("MgnifyClient")
# .mgnify_get_x_for_y(cl, "MGYS00005126", "studies", "samples")
.mgnify_get_x_for_y <- function(client, x, typeX, typeY, use.cache = FALSE, ...){
    # Fetch the data on samples/analyses as a json list
    res <- .mgnify_retrieve_json(
        client,
        paste(typeX, x, sep = "/"),
        use.cache = use.cache,
        ...)
    # Get related analyses when samples were found and vice versa if result was found.
    if( !is.null(res) ){
        res <- res[[1]]$relationships[[typeY]]$links$related
    }
    return(res)
}

############################## .mgnify_get_x_for_y #############################
# Internal function to actually perform the http request. Build up the URL then
# issues a GET, parsing the returned JSON into a nested list (uses jsonlite
# internally?) Previously cached results may be retrieved from disk without
# resorting to calling the MGnify server.

# Low level MGnify API handler
#
# .mgnify_retrieve_json deals with handles the actual HTTP GET calls for the
# MGnifyR package, handling API pagination, local result caching, and
# authentication cookies for access to restricted or pre-release datasets.
# Although principally intended for internal MGnifyR use, it's exported for
# direct invocation. Generally though it's not recommended for use by users.
#
# Parameters:
# client MGnifyR client
# path top level search point for the query. One of biomes, samples, runs etc.
# Basically includes all parts of the URL between the base API url and the
# parameter specifications
# complete_url complete url to search, usually retrieved from previous query in
# the "related" section.
# qopts named list or vector containing options/filters to be URL encoded and
# appended to query as key/value pairs
# max.hits Maximum number of data entries to return. The actual number of hits
# returned may be higher than this value, as this parameter only clamps after
# each full page is processed. Set to <=0 to disable - i.e. retrieve all items.
# use.cache Should successful queries be cached on disk locally? There are
# unresolved questions about whether this is a sensible thing to do, but it
# remains as an option. It probably makes sense for single accession grabs,
# but not for (filtered) queries - which are liable to change as new data is
# added to MGnify. Also caching only works for the first page.
# Debug Should we print out lots of information while doing the grabbing?
#
# Return:
# list of results after pagination is dealt with.

#' @importFrom urltools parameters parameters<-
#' @importFrom httr add_headers
#' @importFrom httr GET
#' @importFrom httr config
#' @importFrom httr content
.mgnify_retrieve_json <- function(
        client, path = "biomes", complete_url = NULL, qopts = NULL,
        max.hits = 200, use.cache = FALSE, Debug=FALSE, ...){
    # Warning message if data is not found
    warning_msg <- paste0(path, ": No data found.")
    # client@warnings turns on debugging too:
    if(client@warnings){
        Debug <- TRUE
    }
    # Set up the base url
    # Are we using internal paths?
    if (is.null(complete_url)){
        fullurl <- paste(client@url, path, sep="/")
    } else{
        # Or direct links from e.g. a "related" section
        # Set the full url, but clean off any existing parameters
        # (page, format etc) as they'll be added back later:
        fullurl <- complete_url
        parameters(fullurl) <- NULL
        path <- substr(fullurl, nchar(client@url) + 2, nchar(fullurl))
    }

    # Convert to csv if filters are lists.
    # This doesn't check if they  can  be searched for in the API,
    # which is an issue since no error is returned by the JSON if the search
    # is invalid - we only get a result as if no query was present...
    tmpqopts <- lapply(qopts,function(x) paste(x,collapse = ','))

    # Include the json and page position options
    # full_qopts <- as.list(c(format="json", tmpqopts, page=1))
    full_qopts <- as.list(c(format="json", tmpqopts))

    # Build up the cache name anyway - even if it's not ultimately used:
    fname_list <- c(path, names(unlist(full_qopts)), unlist(full_qopts))
    cache_fname <- paste(fname_list,collapse = "_")
    cache_full_fname <- paste(client@cacheDir, '/', cache_fname, '.RDS', sep="")

    # Quick check to see if we should clear the disk cache  for this
    # specific call  - used for debugging and when MGnify breaks
    if(use.cache & client@clearCache){
        if( file.exists(cache_full_fname) ){
            message(paste("clear_cache is TRUE: deleting ", cache_full_fname, sep=""))
            unlink(cache_full_fname)
        }
    }

    # Do we want to try and use a cache to speed things up?
    if(use.cache & file.exists(cache_full_fname)){
        final_data <- readRDS(cache_full_fname)
    } else{
        # Authorization: Bearer <your_token>
        if(!is.null(client@authTok)){
            add_headers(
                .headers = c(Authorization = paste("Bearer",
                                                   client@authTok, sep=" ")))
        }
        res <- GET(url=fullurl, config(verbose=Debug), query=full_qopts )
        # Get the data
        data <- content(res, ...)

        # Check if the search was successful and data can be found
        not_found <- (res$status_code != 200) || (
            is.null(data$data) || length(data$data) == 0)
        # If data is found
        if( !not_found ){
            # Fetch all the data
            final_data <- .retrieve_json_data(
                client, data, fullurl, full_qopts, max.hits, Debug
            )
        } else{
            final_data <- NULL
            if( res$status_code != 200 ){
                warning_msg <- paste0(path, ": ", data$errors[[1]]$detail)
            }
        }
        # Save the result to file if specified
        if (use.cache && !file.exists(cache_full_fname)){
            # Make sure the directory is created...
            dir.create(dirname(cache_full_fname), recursive = TRUE,
                       showWarnings = client@warnings)
            saveRDS(final_data, file = cache_full_fname)
        }
    }
    # Give warning if data is not found.
    if( is.null(final_data) ){
        warning(warning_msg, call. = FALSE)
    }
    return(final_data)
}

# This retrives all the data related to accession. FOr example, it loops
# oer multiple pages.
.retrieve_json_data <- function(
        client, data, fullurl, full_qopts, max.hits, Debug, ...){
    # At this point, data$data is either a list of lists or a single named
    # list. If it's a single entry, it needs embedding in a list for
    # consistency downstream datlist is built up as a list of pages, where
    # each entry must be another list. Thus, on the first page,
    #
    datlist <- list()
    if( !is.null(names(data$data)) ){
        # Create something to store the returned data
        datlist[[1]] <- list(data$data)
    }else{
        datlist[[1]] <- data$data
    }
    # Check to see if there's pagination required
    if( "meta" %in% names(data) ){
        # Yes, paginate
        pstart <- as.numeric(data$meta$pagination$page)
        pend <- as.numeric(data$meta$pagination$pages)
        # We've already got the first one
        if( pend > 1 ){
            # Loop over pages and save their result to list
            for (p in seq(pstart+1,pend)){
                full_qopts$page <- p
                if(!is.null(client@authTok)){
                    add_headers(
                        .headers = c(
                            Authorization = paste("Bearer", client@authTok, sep=" ")))
                }
                curd <- content(GET(fullurl, config(verbose=Debug),
                                    query=full_qopts ), ...)
                datlist[[p]] <- curd$data
                # Check to see if we've pulled enough entries.
                # With NULL and -1, disable max.hits
                curlen <- sum(sapply(datlist, length))
                if( !is.null(max.hits) && curlen >= max.hits && max.hits != -1 ){
                    break
                }
            }
        }
    }
    # Combine results from different pages
    final_data <- unlist(datlist, recursive = FALSE)
    return(final_data)
}

#Internal functions to parse the attributes/hierarchy list into a data.frame
.mgnify_parse_tax <- function(json){
    df <- as.data.frame(
        c(json$attributes["count"], unlist(json$attributes$hierarchy)),
        stringsAsFactors = FALSE)
    df$index_id <- json$attributes$lineage
    df

}
.mgnify_parse_func <- function(json){
    df <- as.data.frame(json$attributes, stringsAsFactors = FALSE)
    df$index_id <- json$attributes$accession
    df
}
