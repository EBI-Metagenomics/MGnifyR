#'Download arbitray files from MGnify, including processed reads and identified protein sequences.
#'
#'\code{mgnify_download} is a convenient wrapper round generic the url downloading functionality in R, taking care of things like local
#'caching and authentication. By default, \code{mgnify_download}
#'
#' @importFrom httr add_headers
#' @importFrom httr content
#' @importFrom httr write_disk
#'
#'@param client MGnify client object
#'@param url The url of the file we wish to download
#'@param target_filename An optional local filename to use for saving the file. If NULL (default), MGNify local cache settings will be used.
#'If the file is intended to be processed in a seperate program, it may be sensible to provide a meaningful \code{target_filename}, rather than having to hunt
#'through the cache folders. If \code{target_filename} is NULL \emph{and} \code{usecache} is \code{FALSE}, the \code{read_func} parameter must be supplied or the file
#'will be downloaded and then deleted.
#'@param read_func An optional function name to process the downloaded file and return the results, rather than relying on post processing. The primary use=case for
#' this parameter is when local disk space is limited and downloaded files can be quickly processed and discarded. The function should take a single parameter,
#' the downloaded filename, and may return any valid R object.
#'@param usecache whether to enable the default MGnifyR caching mechanism. File locations are overridden if \code{target_filename} is supplied.
#'@param Debug whether to enable debug output of the HTTP call - only useful for development.
#'@return Either the local filename of the downloaded file, be it either the location in the MGNifyR cache, or target_filename. If \code{read_func} is used, its result
#' will be returned.
#'@examples
#' #Make a client object
#' mg <- mgnify_client(cache_dir="/tmp/mgcache")
#' #create a vector of accession ids - these happen to be \code{analysis} accessions
#' accession_vect <- c("MGYA00563876", "MGYA00563877", "MGYA00563878", "MGYA00563879", "MGYA00563880" )
#' downloads <- mgnify_get_downloads(mg, accession_vect, "analyses")
#'
#' #Filter to find the urls of 16S encoding sequences
#' url_list <- downloads[downloads$attributes.description.label == "Contigs encoding SSU rRNA","download_url"]
#'
#' #Example 1:
#' #Download the first file
#' supplied_filename = mgnify_download(mg, url_list[[1]], target_filename="SSU_file.fasta.gz")
#'
#'
#' #Example 2:
#' #Just use local caching
#' cached_filename = mgnify_download(mg, url_list[[2]])
#'
#' #Example 3:
#' #Using read_func to open the reads with readDNAStringSet from \code{biostrings}. Without retaining on disk
#' dna_seqs <- mgnify_download(mg, url_list[[3]], read_func=readDNAStringSet, usecache=F)
#'
#' @export
mgnify_download <- function(client, url, target_filename=NULL, read_func=NULL, usecache=TRUE, Debug=FALSE){
    #Set up filenames for storing the data
    ftgt=NULL
    if (! is.null(target_filename)){
        file_tgt = target_filename
    }else if(usecache == TRUE){
        #Build a filename out of the url, including the full paths. Annoying, but some downloads (e.g. genome results) are just names like
        # core_genes.fa , which would break the caching.
        cachetgt = gsub(paste(client@url,'/',sep=""), '', url)
        #Make sure the direcory exists

        cache_full_name = paste(client@cache_dir, cachetgt, sep="/")
        dir.create(dirname(cache_full_name), recursive = T, showWarnings = client@warnings)


        file_tgt = cache_full_name
    } else{
        file_tgt = tempfile()[[1]]
    }

    if(usecache & client@clear_cache){
        message(paste("clear_cache is TRUE: deleting ", file_tgt, sep=""))
        tryCatch(unlink(file_tgt), error=warning)
    }

    #Only get the data if it's not already on disk
    if(!(usecache & file.exists(file_tgt))){

        if(!is.null(client@authtok)){
            httr::add_headers(.headers = c(Authorization = paste("Bearer", client@authtok, sep=" ")))
        }
        #If there's an error we need to make sure the cache file isn't written - by default it seems it is.
        tryCatch(expr = {
            curd = httr::content(httr::GET(url, httr::write_disk(file_tgt, overwrite = T)))
        }, error=function(x){
            unlink(file_tgt)
            message(paste("Error retrieving file",file_tgt))
            message(paste("Error:",x))
            stop()
        })
    }

    if (is.null(read_func)){
        result = file_tgt
    } else{
        result = read_func(file_tgt)
    }

    if (is.null(target_filename) & !usecache){
        #Need to clear out the temporary file
        unlink(file_tgt)
    }
    result
}
