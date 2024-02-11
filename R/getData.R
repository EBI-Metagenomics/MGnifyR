#' Versatile function to retrieve results from MGnify database
#'
#' @details
#'
#' @param x A \code{MgnifyClient} object.
#'
#' @param accession A single character value or a vector of character values
#' specifying accession IDs to return results for.
#'
#' @param ... optional arguments:
#'
#' @return
#'
#' @examples
#' # Create a client object
#' mg <- MgnifyClient(useCache = FALSE)
#'
#' @name getData
NULL

#' @rdname getData
#' @include allClasses.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "getData", signature = c(x = "MgnifyClient"), function(
    x, type, accession.type = NULL, accession = NULL, as.df = TRUE, ...){
    ############################### INPUT CHECK ################################
    if( !(.is_non_empty_character(accession)) || is.null(accession) ){
        stop(
            "'accession' must be a single character value or list of ",
            "character values specifying the MGnify accession identifier.",
            call. = FALSE)
    }
    if( !(.is_non_empty_character(accession.type)) || is.null(accession.type) ){
        stop(
            "'accession' must be a single character value or list of ",
            "character values specifying the MGnify accession identifier.",
            call. = FALSE)
    }
    ############################# INPUT CHECK END ##############################
    result <- .get_results_as_json_list(x, type, accession.type, accession, ...)
    # Create df
    if( as.df ){
        result <- .convert_json_list_to_df(result)
    } else if( length(result) == 1 ){
        result <- result[[1]]
    }
    return(result)
})

################################ HELP FUNCTIONS ################################

#' @importFrom plyr llply
.get_results_as_json_list <- function(x, type, accession.type, accession, ...){
    # Create a path
    if( !is.null(accession.type) && !is.null(accession) ){
        path = paste0(accession.type, "/", accession, "/", type)
    } else{
        path = type
    }
    # Find results
    res <- llply(path, function(x){
        .mgnify_retrieve_json(mg, path = x, ...)
    })
    return(res)
}

#' @importFrom tidyjson spread_all
#' @importFrom dplyr bind_rows
.convert_json_list_to_df <- function(result){
    # Create data.frames from individual search results
    result <- lapply(result, function(x) as.data.frame(spread_all(x)))
    # Merge individual data.frames to one
    result <- bind_rows(result)
    # Remove duplicate rows
    result <- result[ !duplicated(result), ]
    return(result)
}

# 1. getJSON function --> get accession.type, accession, type
# accession.type --> type of accession ID, can e also NULL, defaukt --> must be specified if accession is specified.
# accession --> accession IDs, can be also NULL (default)
# type --> type of data to search, must nbe specified
# as.df --> return df
# DO not check datatypes,
# getJSON(mg, type, accession.type = NULL, accession = NULL, ...)
