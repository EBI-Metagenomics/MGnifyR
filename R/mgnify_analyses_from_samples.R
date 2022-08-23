#' Look up analysis accession IDs for one or more sample accessions
#'
#' \code{mgnify_analyses_from_samples} Retrieve Analysis accession IDs associated with the supplied Sample accession
#'
#' Helper function to get all analyses associated with the given samples.
#'
#' @importFrom plyr llply
#'
#' @param client \code{mgnify_client} instance
#' @param accession Single sample accession id, or vector/list of accessions for which to retrieve Analyses ids
#' @param usecache Flag to determine whether to re-use/store data on disk, rather than query the server.
#' @return vector of associated Analysis accession ids
#' @examples
#' #Retrieve all analysis ids from samples
#' result <- mgnify_analyses_from_samples(myclient, c("SRS4392730", "SRS4392743"))
#'
#' @export
mgnify_analyses_from_samples <- function(client, accession, usecache=T){
    #analyses_accessions <- sapply(as.list(accession), function(x){
    analyses_accessions <- llply(as.list(accession), function(x){
        accurl <- mgnify_get_x_for_y(client, x, "samples","analyses", usecache = usecache )
        #For some reason, it appears you "sometimes" have to go from study to runs to analyses. Need
        #to query this with the API people...
        if(is.null(accurl)){
            runurl <- mgnify_get_x_for_y(client, x, "samples","runs", usecache = usecache )
            jsondat <- mgnify_retrieve_json(client, complete_url = runurl, usecache = usecache)
            run_accs <- lapply(jsondat, function(y) y$id)
            a_access <- sapply(as.list(run_accs), function(z){
                accurl <- mgnify_get_x_for_y(client, z, "runs","analyses", usecache = usecache )
                jsondat <- mgnify_retrieve_json(client, complete_url = accurl, usecache = usecache)
                # Now... if jsondat is empty, it means we couldn't find an analysis for this run. This is known to occur when
                # an assembly has been harvested (or something like that). There may be other cases as well. Anyway, what we'll do is
                # go try and look for an assembly->analysis entry instead.
                if(length(jsondat) == 0){
                    assemurl <- mgnify_get_x_for_y(client, z, "runs","assemblies", usecache = usecache )
                    jsondat <- mgnify_retrieve_json(client, complete_url = assemurl, usecache = usecache)
                    assemids <- lapply(jsondat, function(x) x$id)
                    if(length(assemids) >0){
                        #Assumes that there's only one assembly ID per run... I hope that's okay.
                        accurl <- mgnify_get_x_for_y(client, assemids[[1]], "assemblies","analyses", usecache = usecache )
                        jsondat <- mgnify_retrieve_json(client, complete_url = accurl, usecache = usecache)
                    }else{
                        #If we've got to this point, I give up - jsut return an empty list...
                        warning(paste("Failed to find an analysis for sample", accession))
                    }
                }
                lapply(jsondat, function(x) x$id)
            })
            unlist(a_access)
        }else{
            jsondat <- mgnify_retrieve_json(client, complete_url = accurl, usecache = usecache)
            #Just need the accession ID
            lapply(jsondat, function(x) x$id)
        }}, .progress="text")
    unlist(analyses_accessions)
}
