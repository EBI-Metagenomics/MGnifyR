#' Download any MGnify files, also including processed reads and
#' identified protein sequences
#'
#' @details
#' \code{getFile} is a convenient wrapper round generic the URL
#' downloading functionality in R, taking care of things like local
#' caching and authentication.
#'
#' @param x A \code{MgnifyClient} object.
#'
#' @param url A single character value specifying the url address of the file
#' we wish to download.
#'
#' @param file A single character value or NULL specifying an
#' optional local filename to use for saving the file. If NULL (default),
#' MGNify local cache settings will be used. If the file is intended to be
#' processed in a separate program, it may be sensible to provide a
#' meaningful \code{file}, rather than having to hunt through the
#' cache folders. If \code{file} is NULL \emph{and} \code{useCache(client)}
#' is \code{FALSE}, the \code{read.func} parameter must be supplied or the
#' file will be downloaded and then deleted.
#' (By default: \code{file = NULL})
#'
#' @param read.func A function specifying an optional function to process the
#' downloaded file and return the results, rather than relying on post
#' processing. The primary use-case for this parameter is when local disk
#' space is limited and downloaded files can be quickly processed and
#' discarded. The function should take a single parameter, the downloaded
#' filename, and may return any valid R object.
#' (By default: \code{read.func = NULL})
#'
#' @param ... Additional arguments; not used currently.
#'
#' @return Either the local filename of the downloaded file, be it either the
#' location in the MGNifyR cache or file. If \code{read.func} is
#' used, its result will be returned.
#'
#' @examples
#' # Make a client object
#' mg <- MgnifyClient(cache_dir="/tmp/mgcache")
#' # Create a vector of accession ids - these happen to be \code{analysis}
#' # accessions
#' accession_vect <- c(
#'     "MGYA00563876", "MGYA00563877", "MGYA00563878", "MGYA00563879",
#'     "MGYA00563880")
#' downloads <- searchFile(mg, accession_vect, "analyses")
#'
#' # Filter to find the urls of 16S encoding sequences
#' url_list <- downloads[
#'     downloads$attributes.description.label == "Contigs encoding SSU rRNA",
#'     "download_url"]
#'
#' # Example 1:
#' # Download the first file
#' supplied_filename <- getFile(
#'     mg, url_list[[1]], file="SSU_file.fasta.gz")
#'
#'
#' # Example 2:
#' # Just use local caching
#' cached_filename <- getFile(mg, url_list[[2]])
#'
#' # Example 3:
#' # Using read.func to open the reads with readDNAStringSet from
#' # \code{biostrings}. Without retaining on disk
#' dna_seqs <- getFile(
#'     mg, url_list[[3]], read.func = readDNAStringSet)
#'
#' @name getFile
NULL

#' @rdname getFile
#' @importFrom httr add_headers content write_disk
#' @include allClasses.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod("getFile", signature = c(x = "MgnifyClient"), function(
        x, url, file = NULL, read.func = NULL, ...){
    ############################### INPUT CHECK ################################
    if( !.is_non_empty_string(url) ){
        stop(
            "'url' must be a single character value specifying ",
            "the url of the file.", call. = FALSE)
    }
    if( !(.is_non_empty_string(file) || is.null(file)) ){
        stop(
            "'file' must be NULL or a single character value ",
            "specifying the name of file being saved.", call. = FALSE)
    }
    if( !(is.function(read.func) || is.null(read.func)) ){
        stop(
            "'read.func' must be a function that is used to process the file ",
            "or NULL.", call. = FALSE)
    }
    ############################# INPUT CHECK END ##############################
    # Get file
    result <- .mgnify_download(
        client = x, url = url, file = file,
        read.func = read.func, ...)
    return(result)
})

#' Listing files available for download
#'
#' @details
#' THe function is a wrapper function allowing easy enumeration
#' of downloads available for a given accession (or list thereof). Returns a
#' single data.frame containing all available downloads and associated
#' metadata, including the url location and description. This can then be
#' filtered to extract the urls of interest, before actually
#' retrieving the files using \code{mgnify_download}
#'
#' @param accession A single character value or a vector of character values
#' specifying accession IDs to return results for.
#'
#' @param type A single character value specifying the type of objects to
#' query. Must be one of the following options: \code{analysis}, \code{samples},
#' \code{studies}, \code{assembly}, \code{genome} or \code{run}.
#' (By default: \code{type = "samples"})
#'
#' @return \code{data.frame} containing all discovered downloads. If
#' multiple \code{accessions} are queried, the \code{accessions} column
#' may to filter the results - since rownames are not set (and wouldn;'t
#' make sense as each query will return multiple items)
#'
#' @examples
#' # Make a client object
#' mg <- MgnifyClient(cache_dir="/tmp/mgcache", useCache = TRUE)
#' # Create a vector of accession ids - these happen to be \code{analysis}
#' # accessions
#' accession_vect <- c(
#'     "MGYA00563876", "MGYA00563877", "MGYA00563878",
#'     "MGYA00563879", "MGYA00563880" )
#' downloads <- mgnify_get_download_urls(mg, accession_vect, "analyses")
#'
#' @name getFile
NULL

#' @rdname getFile
#' @importFrom plyr llply rbind.fill
#' @importFrom urltools parameters parameters<-
#' @include allClasses.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod("searchFile", signature = c(x = "MgnifyClient"), function(
        x, accession, type = c(
            "studies", "samples", "analyses", "assemblies", "genomes", "run"),
        ...
        ){
    ############################### INPUT CHECK ################################
    if( !.is_non_empty_character(accession) ){
        stop(
            "'accession' must be a list of character values specifying ",
            "the MGnify accession identifiers.",
            call. = FALSE)
    }
    if( !(.is_non_empty_string(type)) ){
        stop(
            "'type' must be a single character value specifying ",
            "the type of instance to query.", call. = FALSE)
    }
    type <- match.arg(type, several.ok = FALSE)
    ############################# INPUT CHECK END ##############################
    # Get file urls
    result <- .mgnify_get_download_urls(
        client = x, accession = accession, type = type, ...)
    return(result)
})

################################ HELP FUNCTIONS ################################

# Download the specified files from the database
.mgnify_download <- function(
        client, url, file, read.func, use.cache = useCache(client),
        url.address = databaseUrl(client), cache.dir = cacheDir(client),
        show.warnings = showWarnings(client), clear.cache = clearCache(client),
        auth.tok = authTok(client), ...){
    # Input check
    if( !.is_non_empty_string(url.address) ){
        stop(
            "'url.address' must be a string.", call. = FALSE)
    }
    if( !.is_a_bool(use.cache) ){
        stop(
            "'use.cache' must be a single boolean value.", call. = FALSE)
    }
    if( !.is_non_empty_string(cache.dir) ){
        stop(
            "'cache.dir' must be a string.", call. = FALSE)
    }
    if( !.is_a_bool(show.warnings) ){
        stop(
            "'show.warnings' must be a single boolean value.", call. = FALSE)
    }
    if( !.is_a_bool(clear.cache) ){
        stop(
            "'clear.cache' must be a single boolean value.", call. = FALSE)
    }
    if( !(.is_non_empty_string(auth.tok) || is.null(auth.tok)) ){
        stop(
            "'auth.tok' must be a string or NULL.", call. = FALSE)
    }
    #
    # Set up filenames for storing the data
    if ( !is.null(file) ){
        file_path <- file
    }else if(use.cache){
        # Build a filename out of the url, including the full paths. Annoying,
        # but some downloads (e.g. genome results) are just names like
        # core_genes.fa , which would break the caching.
        cachetgt <- gsub(paste(url.address, "/", sep = ""), "", url)

        # Make sure the directory exists
        cache_full_name <- file.path(cache.dir, cachetgt)
        dir.create(
            dirname(cache_full_name), recursive = TRUE,
            showWarnings = show.warnings)
        file_path <- cache_full_name
    } else{
        file_path <- tempfile()[[1]]
    }

    # Clear cache if specified
    if( use.cache && clear.cache && file.exists(file_path) ){
        message("clear_cache is TRUE: deleting ", file_path)
        unlink(file_path)
    }

    # Only get the data if it's not already on disk
    if( !file.exists(file_path) || (use.cache && file.exists(file_path)) ){
        # Add authentication details to query options
        if( !is.null(auth.tok) ){
            add_headers(.headers = c(
                Authorization = paste("Bearer", auth.tok, sep = " ")))
        }
        # If there's an error we need to make sure the cache file isn't written
        # - by default it seems it is.
        res <- GET(url, write_disk(file_path, overwrite = TRUE))
        # If the file was not successfully downloaded
        if( res$status_code != 200 ){
            # Remove the downloaded file
            unlink(file_path)
            stop(
                url, ": ", content(res, ...)$errors[[1]]$detail,
                " Could not load the file from database.",
                call. = FALSE)
        }
    }
    # Whether to use user-specified read function
    if( is.null(read.func) ){
        result <- file_path
    } else{
        result <- read.func(file_path)
    }
    return(result)
}

# Get URL addresses of downloadable files that are related to certain accession
# ID.
.mgnify_get_download_urls <- function(
        client, accession, type, use.cache = useCache(client),
        show.messages = verbose(client), ...){
    # Input check
    if( !.is_a_bool(use.cache) ){
        stop(
            "'use.cache' must be a single boolean value.", call. = FALSE)
    }
    if( !.is_a_bool(show.messages) ){
        stop(
            "'show.messages' must be a single boolean value.", call. = FALSE)
    }
    show.messages <- ifelse(show.messages, "text", "none")
    #
    # Give message about progress
    if( show.messages == "text" ){
        message("Searching files...")
    }
    # L
    # Loop though accession IDs and find the info
    results <- llply(accession, function(x){
        # Get the data as nested json list
        download_list <- .mgnify_retrieve_json(
            client, paste(type, x, "downloads", sep = "/"),
            use.cache = use.cache, ...)
        # Convert to df
        df <- do.call(rbind.fill, lapply(download_list, function(x){
            as.data.frame(x, stringsAsFactors = FALSE)}
            ))
        # Add info to df
        df$accession <- x
        df$type <- type
        # If no match, df is a list --> convert to data.frame
        if( !is.data.frame(df) ){
            df <- as.data.frame(df)
        } else {
            # If search result was found, modify
            # For convenience, rename the "self" column to "download_url" -
            # which is what it actually is...
            colnames(df)[colnames(df) == "self"] <- "download_url"
            # Finally, strip off any options from the url - they sometimes seem
            # to get format=json stuck on the end
            urls <- df$download_url
            parameters(urls) <- NULL
            df$download_url <- urls
        }
        return(df)
    }, .progress = show.messages)
    # Combine results of multiple accessions IDs
    results <- do.call(rbind.fill, results)
    return(results)
}
