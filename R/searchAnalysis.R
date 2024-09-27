#' Look up analysis accession IDs for one or more study or sample accessions
#'
#' @details
#' Retrieve analysis accession IDs associated with the supplied study or
#' sample accession.  In MGnify, an analysis accession refers to a certain
#' pipeline analysis, such as specific 16S rRNA or shotgun metagenomic mapping.
#' Studies can include multiple samples, and each sample can undergo multiple
#' analyses using these pipelines. Each analysis is identified by a unique
#' accession ID, allowing precise tracking and retrieval of analysis results
#' within the MGnify database.
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
#' @return Vector of analysis accession IDs.
#'
#' @examples
#' # Create a client object
#' mg <- MgnifyClient(useCache = FALSE)
#'
#' # Retrieve analysis ids from study MGYS00005058
#' result <- searchAnalysis(mg, "studies", c("MGYS00005058"))
#'
#' \dontrun{
#' # Retrieve all analysis ids from samples
#' result <- searchAnalysis(
#'     mg, "samples", c("SRS4392730", "SRS4392743"))
#' }
#'
#' @name searchAnalysis
NULL

#' @rdname searchAnalysis
#' @importFrom plyr llply
#' @include AllClasses.R AllGenerics.R MgnifyClient.R utils.R
#' @export
setMethod("searchAnalysis", signature = c(x = "MgnifyClient"), function(
        x, type, accession, ...){
    ############################### INPUT CHECK ################################
    if( !(length(type) == 1 && type %in% c("samples", "studies")) ){
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
    result <- .mgnify_analyses_from_studies_and_samples(
        client = x, accession = accession, type = type, ...)
    return(result)
})

################################ HELP FUNCTIONS ################################
# Get analysis accessions based on studies or samples. The result is a vector
# of analyses IDs.
.mgnify_analyses_from_studies_and_samples <- function(
        client, accession, type, show.messages = verbose(client), ...){
    # Input check
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
    # Search analyses IDs
    analysis_ids <- .get_all_analyses_ids(
        client, accession, type, "analyses", show.messages = show.messages, ...)
    # Check which study/sample ID resulted to found analysis ID
    not_found <- accession[ !accession %in% names(analysis_ids) ]
    # If user is searching analyses based on samples, we can still try another
    # approach. Sometimes, those "sample" IDs refer to runs instead.
    if( length(not_found) > 0 && type == "samples" ){
        # Finds runs based on samples
        temp <- .get_all_analyses_ids(
            client, accession, "samples", "runs",
            show.messages = show.messages, ...)
        # Create a data.frame that holds all the IDs to book keep matches
        # between IDs.
        id_df <- data.frame(sample = names(temp), run = temp)
        # Based on those runs, search analyses
        temp <- .get_all_analyses_ids(
            client, id_df[["run"]], "runs", "analyses",
            show.messages = show.messages, ...)
        # Add found analysis IDs to data.frame
        temp_df <- id_df[match(names(temp), id_df[["run"]]), ]
        temp_df[["analyses"]] <- temp
        id_df <- merge(id_df, temp_df, all = TRUE)
        
        # If there still are samples that were not found, we can try to get
        # analyses from assemblies. That is why we try to first fetch assemblies
        # based on runs.
        temp <- .get_all_analyses_ids(
            client, id_df[is.na(id_df[["analyses"]]), "run"], "runs",
            "assemblies", show.messages = show.messages, ...)
        # Add found analysis IDs to data.frame
        temp_df <- id_df[match(names(temp), id_df[["run"]]), ]
        temp_df[["assemblies"]] <- temp
        id_df <- merge(id_df, temp_df, all = TRUE)
        # Then based on assemblies, we can finally try to find analyses.
        temp <- .get_all_analyses_ids(
            client, id_df[is.na(id_df[["analyses"]]), "assemblies"],
            "assemblies", "analyses", show.messages = show.messages, ...)
        # Add found analysis IDs to data.frame
        temp_df <- id_df[match(names(temp), id_df[["assemblies"]]), ]
        temp_df[["analyses"]] <- temp
        id_df <- merge(id_df, temp_df, all = TRUE)
        # Now we should have a table that contains all the analyses that were
        # possible to find. Add these analyses to the original result list.
        temp <- id_df[["analyses"]]
        names(temp) <- id_df[["sample"]]
        temp <- temp[ !is.na(temp) ]
        analysis_ids <- c(analysis_ids, temp)
        # Update the "not found samples" vector
        not_found <- accession[ !accession %in% names(analysis_ids) ]
    }
    # If the data was not found for specified ID, give warning
    if( length(not_found) > 0 ){
        warning(
            "\nAnalyses not found for the following ", type, ": '",
            paste(not_found, collapse = "', '"), "'", call. = FALSE)
    }
    return(analysis_ids)
}

# This function gets IDs type "type_from" as input and tries to fetch
# corresponding IDs type "type_to".
# based on those studies or samples.
.get_all_analyses_ids <- function(
        client, ids, type_from, type_to, show.messages,
        use.cache = useCache(client), ...){
    #
    if( !.is_a_bool(use.cache) ){
        stop(
            "'use.cache' must be a single boolean value", call. = FALSE)
    }
    #
    # Get only unique IDs
    ids <- unique(ids)
    # Loop through accessions
    analysis_ids <- llply(ids, function(id){
        # Get URL address of results that were found. For instance, URL address
        # of analyses based on study ID/accession
        url <- .mgnify_get_x_for_y(
            client, id, type_from, type_to, use.cache = use.cache,
            ...)
        # Check whether results were found or not
        res <- NULL
        if( !is.null(url) ){
            # Get data
            json <- .mgnify_retrieve_json(
                client, complete_url = url, use.cache = use.cache,
                max.hits = NULL, ...)
            # We need just the accession ID
            res <- lapply(json, function(x) x$id) |> unlist()
            # Add accession as name. There might be multiple analyses for each
            # accession. This helps to determine which analyses belong to which
            # study.
            if( length(res) > 0 ){
                names(res) <- rep(id, length(res))
            }
        }
        return(res)
    }, .progress = show.messages)
    # Create a vector from results
    analysis_ids <- analysis_ids |> unlist()
    return(analysis_ids)
}
