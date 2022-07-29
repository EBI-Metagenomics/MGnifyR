#' Listing files available for download
#'
#' \code{mgnify_get_downloads} is a wrapper function allowing easy enumeration of downloads available for a given
#' accession (or list thereof). Returns a single data.frame containing all available downloads and associated metadata,
#' including the url location and description. This can then be filtered to extract the urls of interest, before actually
#' retrieving the files using \code{mgnify_download}
#'
#' @importFrom plyr llply
#' @importFrom plyr rbind.fill
#' @importFrom urltools parameters
#'
#'@param client valid MGnify client object
#'@param accessions list of accessions to query
#'@param accession_type one of \code{analysis},\code{samples},\code{studies},\code{assembly},\code{genome} or \code{run}
#'@param usecache whether to use the on-disk cache to speed up queries (default T)
#'@return \code{data.frame} containing all discovered downloads. If multiple \code{accessions} are queried, the \code{accessions} column
#' may to filter the results - since rownames are not set (and wouldn;'t make sense as each query will return multiple items)

#'@examples
#' #Make a client ibject
#' mg <- mgnify_client(cache_dir="/tmp/mgcache")
#' #create a vector of accession ids - these happen to be \code{analysis} accessions
#' accession_vect <- c("MGYA00563876", "MGYA00563877", "MGYA00563878", "MGYA00563879", "MGYA00563880" )
#' downloads <- mgnify_get_downloads(mg, accession_vect, "analyses")
#'@export

mgnify_get_download_urls <- function(client, accessions, accession_type, usecache=T){
    results <- plyr::llply(accessions, function(x){
        download_list <- mgnify_retrieve_json(client, paste(accession_type,x,"downloads", sep="/"), usecache = usecache)
        df <- do.call(plyr::rbind.fill,lapply(download_list, function(x) as.data.frame(x,stringsAsFactors=F)))
        df$accession <- x
        df$accession_type <- accession_type
        #for convenience, rename the "self" column to "download_url" - which is what it actually is...
        colnames(df)[colnames(df) == 'self'] <- 'download_url'
        #finally, strip off any options from the url - they sometimes seem to get format=json stuck on the end
        urls <- df$download_url
        urltools::parameters(urls) <- NULL
        df$download_url <- urls
        df
    }, .progress="text")
    do.call(plyr::rbind.fill, results)
}
