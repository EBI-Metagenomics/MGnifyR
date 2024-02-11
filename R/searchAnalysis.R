#' Look up analysis accession IDs for one or more study or sample accessions
#'
#' @details
#' Retrieve analysis accession IDs associated with the supplied study or
#' sample accession.
#'
#' @param x A \code{MgnifyClient} object.
#'
#' @param type A single character value specifying a type of
#' accession IDs specified by \code{accession}. Must be "studies" or "samples".
#'
#' @param accession A single character value or a vector of character values
#' specifying study or sample accession IDs that are used to retrieve analyses
#' IDs.
#'
#' @param ... Optional arguments; not currently used.
#'
#' @return vector of analysis accession IDs.
#'
#' @examples
#' # Create a client object
#' mg <- MgnifyClient(useCache = FALSE)
#'
#' # Retrieve all analysis ids from studies
#' # MGYS00005058, MGYS00005058 and MGYS00005058
#' result <- searchAnalysis(mg, "studies", c("MGYS00005058"))
#'
#' \donttest{
#' # Retrieve all analysis ids from samples
#' result <- searchAnalysis(
#'     mg, "samples", c("SRS4392730", "SRS4392743"))
#' }
#'
#' @name searchAnalysis
NULL

#' @rdname searchAnalysis
#' @importFrom plyr llply
#' @include allClasses.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod("searchAnalysis", signature = c(x = "MgnifyClient"), function(
        x, type, accession, ...){
    ############################### INPUT CHECK ################################
    if( !(type %in% c("samples", "studies")) ){
        stop(
            "'type' must be 'samples' or 'studies'.", call. = FALSE)
    }
    if( !(.is_non_empty_character(accession)) ){
        stop(
            "'accession' must be a single character value or vector of ",
            "character values specifying the MGnify accession identifier.",
            call. = FALSE)
    }
    ############################# INPUT CHECK END ##############################
    # Get analysis accession IDs based on sample or study accessions
    if( type == "samples" ){
        result <- .mgnify_analyses_from_samples(
            client = x, accession = accession, ...)
    } else{
        result <- .mgnify_analyses_from_studies(
            client = x, accession = accession, ...)
    }
    return(result)
})

################################ HELP FUNCTIONS ################################
# Get analysis accessions based on studies
.mgnify_analyses_from_studies <- function(
        client, accession, use.cache = useCache(client),
        show.messages = verbose(client), ...){
    # Input check
    if( !.is_a_bool(use.cache) ){
        stop(
            "'use.cache' must be a single boolean value", call. = FALSE)
    }
    if( !.is_a_bool(show.messages) ){
        stop(
            "'show.messages' must be a single boolean value.", call. = FALSE)
    }
    show.messages <- ifelse(show.messages, "text", "none")
    #
    # Give message about progress
    if( show.messages == "text" ){
        message("Fetching analyses...")
    }
    # Loop over studies, get analyses accessions
    analyses_accessions <- llply(as.list(accession), function(x){
        # Find analyses based on studies. Get URL address.
        accurl <- .mgnify_get_x_for_y(
            client, x, "studies","analyses", use.cache = use.cache, ...)
        # If found
        if( !is.null(accurl) ){
            # Get data
            jsondat <- .mgnify_retrieve_json(
                client, complete_url = accurl, use.cache = use.cache,
                max.hits = NULL, ...)
            # Just need the accession ID
            res <- lapply(jsondat, function(x) x$id)
        } else {
            res <- accurl
            warning("Analyses not found for studies ", x, call. = FALSE)
        }
        return(res)
    }, .progress=show.messages)
    res <- unlist(analyses_accessions)
    return(res)
}

# Get analysis accessions based on sample accessions
.mgnify_analyses_from_samples <- function(
        client, accession, use.cache = useCache(client),
        show.messages = verbose(client), ...){
    # Input check
    if( !.is_a_bool(use.cache) ){
        stop(
            "'use.cache' must be a single boolean value", call. = FALSE)
    }
    if( !.is_a_bool(show.messages) ){
        stop(
            "'show.messages' must be a single boolean value.", call. = FALSE)
    }
    show.messages <- ifelse(show.messages, "text", "none")
    #
    # Give message about progress
    if( show.messages == "text" ){
        message("Fetching analyses...")
    }
    # Loop over sample accessions
    analyses_accessions <- llply(as.list(accession), function(x){
        accurl <- .mgnify_get_x_for_y(
            client, x, "samples", "analyses", use.cache = use.cache, ...)
        # For some reason, it appears you "sometimes" have to go from study
        # to runs to analyses. Need to query this with the API people...
        if( is.null(accurl) ){
            temp <- .mgnify_analyses_from_samples_based_on_runs(
                client, x, use.cache, ...)
        } else {
            jsondat <- .mgnify_retrieve_json(
                client, complete_url = accurl, use.cache = use.cache, ...)
            # Just need the accession ID
            temp <- lapply(jsondat, function(x) x$id)
        }
        return(temp)
        }, .progress = show.messages)
    res <- unlist(analyses_accessions)
    return(res)
}

# Get analysis accessions based on runs or assemblies
.mgnify_analyses_from_samples_based_on_runs <- function(
        client, x, use.cache = useCache(client), ...){
    # Input check
    if( !.is_a_bool(use.cache) ){
        stop(
            "'use.cache' must be a single boolean value", call. = FALSE)
    }
    #
    # Get urö for runs
    runurl <- .mgnify_get_x_for_y(
        client, x, "samples","runs", use.cache = use.cache, ...)
    if(is.null(runurl)){
        warning("Analyses not found for samples ", x, call. = FALSE)
        return(runurl)
    }
    # If found, get data for runs
    jsondat <- .mgnify_retrieve_json(
        client, complete_url = runurl, use.cache = use.cache, ...)
    # Get accession ID for the runs
    run_accs <- lapply(jsondat, function(y) y$id)
    # Loop through runs
    analyses_accessions <- lapply(as.list(run_accs), function(z){
        # Get data url of related analyses
        accurl <- .mgnify_get_x_for_y(
            client, z, "runs","analyses", use.cache = use.cache, ...)
        # Get data of those analyses
        jsondat <- .mgnify_retrieve_json(
            client, complete_url = accurl, use.cache = use.cache, ...)
        # Now... if jsondat is empty, it means we couldn't find an
        # analysis for this run. This is known to occur when an assembly
        # has been harvested (or something like that). There may be
        # other cases as well. Anyway, what we'll do is go try and look
        # for an assembly->analysis entry instead.
        if(length(jsondat) == 0){
            # Get url addresses for assemblies based on runs
            assemurl <- .mgnify_get_x_for_y(
                client, z, "runs","assemblies", use.cache = use.cache, ...)
            # Get data on those assemblies
            jsondat <- .mgnify_retrieve_json(
                client, complete_url = assemurl, use.cache = use.cache, ...)
            # Get accession IDs for assemblies
            assemids <- lapply(jsondat, function(x) x$id)
            if(length(assemids) > 0){
                # Assumes that there's only one assembly ID per run...
                # I hope that's okay.
                # Get analyses based on assemblies
                accurl <- .mgnify_get_x_for_y(
                    client, assemids[[1]], "assemblies", "analyses",
                    use.cache = use.cache, ...)
                # Get the data on analyses
                jsondat <- .mgnify_retrieve_json(
                    client, complete_url = accurl, use.cache = use.cache, ...)
            } else{
                # If we've got to this point, I give up - just return an empty
                # list...
                warning(
                    "Failed to find an analysis for sample ", x, call. = FALSE)
            }
        }
        # Get analyses IDs
        if( !is.null(jsondat) ){
            temp <- lapply(jsondat, function(x) x$id)
        } else{
            temp <- NULL
        }
        return(temp)
    })
    analyses_accessions <- unlist(analyses_accessions)
}
