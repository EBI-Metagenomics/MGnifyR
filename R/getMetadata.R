#' Get all Study, Sample and Analysis metadata for the supplied analyses
#' accessions
#'
#' @details
#' The function retrieves all associated study, sample and analysis
#' metadata attributes as a list of analyses accessions.
#'
#' @param x A \code{MgnifyClient} object.
#'
#' @param accession A single character value or a vector of analysis accession
#' IDs specifying accessions to retrieve data for.
#' 
#' @param ... Optional arguments; not currently used.
#'
#' @return \code{data.frame} of metadata for each analysis in the
#' \code{accession} list.
#'
#' @examples
#' # Create a client object
#' mg <- MgnifyClient(useCache = FALSE)
#'
#' # Download all associated study/sample and analysis metadata
#' accession_list <- c("MGYA00377505")
#' meta_dataframe <- getMetadata(mg, accession_list)
#'
#' @name getMetadata
NULL

#' @rdname getMetadata
#' @importFrom plyr llply
#' @importFrom dplyr bind_rows
#' @include AllClasses.R AllGenerics.R MgnifyClient.R utils.R
#' @export
setMethod("getMetadata", signature = c(x = "MgnifyClient"), function(
        x, accession, ...){
    ############################### INPUT CHECK ################################
    if( !is.character(accession) ){
        stop(
            "'accession' must be a single character or a list of character ",
            "values.", call. = FALSE)
    }
    ############################# INPUT CHECK END ##############################
    # Get metadata
    result <- .mgnify_get_analyses_metadata(
        client = x, accession = accession, ...)
    return(result)
})

################################ HELP FUNCTIONS ################################

# Fetch metadata based on analysis accessions.
.mgnify_get_analyses_metadata <- function(
        client, accession, use.cache = useCache(client),
        show.messages = verbose(client), ...){
    # Input check
    if( !.is_a_bool(use.cache) ){
        stop(
            "'use.cache' must be a single boolean value specifying whether to ",
            "show progress.", call. = FALSE)
    }
    if( !.is_a_bool(show.messages) ){
        stop(
            "'show.messages' must be a single boolean value.", call. = FALSE)
    }
    show.messages <- ifelse(show.messages, "text", "none")
    #
    # Give message about progress
    if( show.messages == "text" ){
        message("Fetching metadata...")
    }
    # Loop through analysis accessions and find metadata
    reslist <- llply(as.list(accession), function(x){
        .mgnify_get_single_analysis_metadata(
            client, x, use.cache = use.cache, ...)
    }, .progress = show.messages)
    # Combine all metadata to single df
    df <- do.call(bind_rows, reslist)
    return(df)
}

# Retrieves combined study/sample/analysis metadata - not exported
.mgnify_get_single_analysis_metadata <- function(
        client, accession, use.cache = useCache(client), max.hits = NULL, ...){
    # Input check
    if( !.is_a_bool(use.cache) ){
        stop(
            "'use.cache' must be a single boolean value specifying whether to ",
            "show progress.", call. = FALSE)
    }
    #
    # Get data in json format
    dat <- .mgnify_retrieve_json(
        client, paste("analyses", accession, sep="/"), use.cache = use.cache,
        max.hits = max.hits, ...)
    # If metadata was not found, return the NULL value
    if(is.null(dat)){
        warning(
            "\nFailed to find study metadata for ", accession, call. = FALSE)
        return(dat)
    }

    # There  should  be just a single result
    top_data <- dat[[1]]
    # Convert hit result to df
    analysis_df <- .mgnify_attr_list_to_df_row(
        top_data, metadata_key = "analysis-summary")

    # Build up the metadata dataframe from the analyses_metadata_headers vector:
    sample_met <- .mgnify_retrieve_json(
        client, complete_url = top_data$relationships$sample$links$related,
        use.cache = use.cache, ...)
    study_met <- .mgnify_retrieve_json(
        client, complete_url = top_data$relationships$study$links$related,
        use.cache = use.cache, ...)
    # Again, convert to df
    if(!is.null(sample_met)){
        sample_df <- .mgnify_attr_list_to_df_row(
            sample_met[[1]], metadata_key = "sample-metadata")
    } else{
        warning(
            "\nFailed to find sample metadata for ", accession, call. = FALSE)
        sample_df <- data.frame(accession=NA)
    }
    # It turns out that a sample might not be part of a study - if it's been
    # harvested...
    if(!is.null(study_met)){
        study_df <- .mgnify_attr_list_to_df_row(study_met[[1]])
    } else{
        warning(
            "\nFailed to find study metadata for ", accession, call. = FALSE)
        study_df <- data.frame(accession=NA)
    }
    # Add colnames to sample, study and analysis tables
    colnames(sample_df) <- paste("sample", colnames(sample_df), sep="_")
    colnames(study_df) <- paste("study", colnames(study_df), sep="_")
    colnames(analysis_df) <- paste("analysis", colnames(analysis_df), sep="_")
    # Add what analysis corresponds what sample and study
    rownames(sample_df) <- rownames(analysis_df)
    rownames(study_df) <- rownames(analysis_df)
    # Combine sample and study result
    full_df <- cbind(analysis_df, study_df, sample_df)

    # Extras - include some more metadata from various places
    # Assembly accession
    if("id" %in% names(top_data$relationships$assembly$data)){
        full_df$assembly_accession <- top_data$relationships$assembly$data$id
    }
    # Run accession
    if("id" %in% names(top_data$relationships$run$data)){
        full_df$run_accession <- top_data$relationships$run$data$id
    }
    # biom (from the sample metadata)
    if( !is.null(sample_met[[1]]$relationships$biome$data$id) ){
        full_df$biome_string <- sample_met[[1]]$relationships$biome$data$id
    } else {
        warning("\nFailed to find biome entry for ", accession, call = FALSE)
    }
    return(full_df)
}
