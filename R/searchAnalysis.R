#' Look up analysis accession IDs for one or more study or sample accessions
#'
#' @details
#' Retrieve Analysis accession IDs associated with the supplied study or
#' sample accession.
#'
#' @param x A \code{MgnifyClient} object.
#'
#' @param type A single character value specifying a type of
#' accession IDs specified by \code{accession}. Must be "study" or "sample".
#'
#' @param accession A single character value or a vector of character values
#' specifying study or sample accession IDs that are used to retrieve analyses
#' IDs.
#'
#' @param use.cache A single boolean value to specify whether to re-use/store
#' data on disk, rather than query the server.
#' (By default: \code{use.cache = TRUE})
#'
#' @param ... Optional arguments; not currently used.
#'
#' @return vector of analysis accession IDs.
#'
#' @examples
#' \dontrun{
#' # Retrieve all analysis ids from studies
#' # MGYS00005058, MGYS00005058 and MGYS00005058
#' result <- searchAnalysis(myclient, "study", c("MGYS00005058"))
#'
#' # Retrieve all analysis ids from samples
#' result <- searchAnalysis(
#'     myclient, "sample", c("SRS4392730", "SRS4392743"))
#' }
#'
#' @name searchAnalysis
NULL

#' @rdname searchAnalysis
#' @include MgnifyClient.R
#' @importFrom plyr llply
#' @export
setGeneric("searchAnalysis", signature = c("x"), function(
        x, type, accession, use.cache=TRUE,
        ...
        )
    standardGeneric("searchAnalysis"))

#' @rdname searchAnalysis
#' @export
setMethod("searchAnalysis", signature = c(x = "MgnifyClient"), function(
        x, type, accession, use.cache=TRUE,
        ...
        ){
    ############################### INPUT CHECK ################################
    if( !(type %in% c("sample", "study")) ){
        stop("'type' must be 'sample' or 'study'.",
             call. = FALSE)
    }
    if( !(.is_non_empty_character(accession)) ){
        stop("'accession' must be a single character value or list of ",
             "character values specifying the MGnify accession identifier.",
             call. = FALSE)
    }
    if( !.is_a_bool(use.cache) ){
        stop("'use.cache' must be a single boolean value specifying whether to ",
             "use on-disk caching.", call. = FALSE)
    }
    ############################# INPUT CHECK END ##############################
    # Get analysis accession IDs based on sample or study accessions
    if( type == "sample" ){
        result <- .mgnify_analyses_from_samples(
            client = x, accession = accession, use.cache = use.cache)
    } else{
        result <- .mgnify_analyses_from_studies(
            client = x, accession = accession, use.cache = use.cache)
    }
    return(result)
})

################################ HELP FUNCTIONS ################################

.mgnify_analyses_from_studies <- function(client, accession, use.cache){
    # Loop over studies, get analyses accessions
    analyses_accessions <- llply(as.list(accession), function(x){
        accurl <- .mgnify_get_x_for_y(
            client, x, "studies","analyses", use.cache = use.cache )
        if( !is.null(accurl) ){
            jsondat <- .mgnify_retrieve_json(
                client, complete_url = accurl, use.cache = use.cache, max.hits = -1)
            # Just need the accession ID
            res <- lapply(jsondat, function(x) x$id)
        } else {
            res <- accurl
            warning("Analyses not found for studies ", x, call. = FALSE)
        }
        return(res)
    }, .progress="text")
    res <- unlist(analyses_accessions)
    return(res)
}

.mgnify_analyses_from_samples <- function(client, accession, use.cache){
    # Loop over sample accessions
    analyses_accessions <- llply(as.list(accession), function(x){
        accurl <- .mgnify_get_x_for_y(
            client, x, "samples", "analyses", use.cache = use.cache )
        # For some reason, it appears you "sometimes" have to go from study
        # to runs to analyses. Need to query this with the API people...
        if(is.null(accurl)){
            runurl <- .mgnify_get_x_for_y(
                client, x, "samples","runs", use.cache = use.cache )
            if(is.null(runurl)){
                warning("Analyses not found for samples ", x, call. = FALSE)
                return(runurl)
            }
            jsondat <- .mgnify_retrieve_json(
                client, complete_url = runurl, use.cache = use.cache)
            run_accs <- lapply(jsondat, function(y) y$id)
            a_access <- sapply(as.list(run_accs), function(z){
                accurl <- .mgnify_get_x_for_y(
                    client, z, "runs","analyses", use.cache = use.cache )
                jsondat <- .mgnify_retrieve_json(
                    client, complete_url = accurl, use.cache = use.cache)
                # Now... if jsondat is empty, it means we couldn't find an
                # analysis for this run. This is known to occur when an assembly
                # has been harvested (or something like that). There may be
                # other cases as well. Anyway, what we'll do is go try and look
                # for an assembly->analysis entry instead.
                if(length(jsondat) == 0){
                    assemurl <- .mgnify_get_x_for_y(
                        client, z, "runs","assemblies", use.cache = use.cache )
                    jsondat <- .mgnify_retrieve_json(
                        client, complete_url = assemurl, use.cache = use.cache)
                    assemids <- lapply(jsondat, function(x) x$id)
                    if(length(assemids) >0){
                        #Assumes that there's only one assembly ID per run...
                        # I hope that's okay.
                        accurl <- .mgnify_get_x_for_y(
                            client, assemids[[1]], "assemblies", "analyses",
                            use.cache = use.cache )
                        jsondat <- .mgnify_retrieve_json(
                            client, complete_url = accurl, use.cache = use.cache)
                    }else{
                        # If we've got to this point, I give up - just return an empty list...
                        warning(paste(
                            "Failed to find an analysis for sample ", accession))
                    }
                }
                lapply(jsondat, function(x) x$id)
            })
            unlist(a_access)
        }else{
            jsondat <- .mgnify_retrieve_json(
                client, complete_url = accurl, use.cache = use.cache)
            # Just need the accession ID
            lapply(jsondat, function(x) x$id)
        }}, .progress="text")
    res <- unlist(analyses_accessions)
    return(res)
}
