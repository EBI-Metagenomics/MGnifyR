#' Look up analysis accession IDs for one or more study accessions
#'
#' \code{mgnify_analyses_from_studies} Retrieve Analysis accession IDs associated with the supplied Study accession
#'
#' Helper function to get all analyses associated with the given studies.
#'
#' @importFrom plyr llply
#'
#' @param client \code{mgnify_client} instance
#' @param accession Single study accession id, or vector/list of accessions for which to retrieve Analyses ids
#' @param usecache Flag to determine whether to re-use/store data on disk, rather than query the server.
#' @return vector of Analysis accession ids
#' @examples
#' #Retrieve all analysis ids from studies MGYS00005058, MGYS00005058 and MGYS00005058
#' result <- mgnify_analyses_from_studies(myclient, c("MGYS00005058"))
#'
#' @export
mgnify_analyses_from_studies <- function(client, accession, usecache=T){
    analyses_accessions <- llply(as.list(accession), function(x){
        accurl <- mgnify_get_x_for_y(client, x, "studies","analyses", usecache = usecache )
        jsondat <- mgnify_retrieve_json(client, complete_url = accurl, usecache = usecache, maxhits = -1)
        #Just need the accession ID
        lapply(jsondat, function(x) x$id)
    }, .progress="text")
    unlist(analyses_accessions)
}
