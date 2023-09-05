#' Search MGnify database for studies, samples and runs
#'
#' @details
#' \code{doQuery} is a flexible query function, harnessing the "full"
#' power of the JSONAPI MGnify search filters. Search results may be filtered
#' by metadata value, associated study/sample/analyese etc. Details of the
#' capabilities may be found
#' [here](https://emg-docs.readthedocs.io/en/latest/api.html#customising-queries).
#' Currently, the following filters are available (based on examination of the
#' Python source code):
#'     \itemize{
#'         \item{\strong{Studies}: accession, biome_name, lineage, centre_name}
#'         \item{\strong{Samples}: accession, experiment_type, biome_name,
#'             lineage, geo_loc_name, latitude_gte, latitude_lte,
#'             longitude_gte, longitude_lte, species, instrument_model,
#'             instrument_platform, metadata_key, metadata_value_gte,
#'             metadata_value_lte, metadata_value, environment_material,
#'             environment_feature, study_accession}
#'        \item{\strong{Runs}: accession, experiment_type, biome_name, lineage,
#'            species, instrument_platform, instrument_model}
#'        }
#' Unfortunately it appears that in some cases, some of these filters don't work
#' as expected, so it is important to check the results returned match up with
#' what's expected. Even more unfortunately if there's an error in the parameter
#' specification, the query will run as if no filter parameters were present
#' at all. Thus the result will appear superficially correct but will infact
#' correspond to something completely different. This behaviour will hopefully
#' be fixed in future incarnations of the MGnifyR or JSONAPI, but for now users
#' should double check returned values.
#'
#' It is currently not possible to combine queries of the same type in a single
#' call (for example to search for samples \emph{between} latitude). However,
#' it is possible to run multiple queries and combine the results using set
#' operations in R to get the desired behaviour.
#'
#' @param x A \code{MgnifyClient} object.
#'
#' @param type A single character value specifying the type of objects to
#' query. Must be one of the following options: \code{studies}, \code{samples},
#' \code{runs} or \code{analyses}. (By default: \code{type = "studies"})
#'
#' @param accession A single character value or a vector of character values
#' specifying MGnify accession identifiers (of type \code{type}) or NULL. When
#' NULL, all results defined by other parameters are retrieved.
#' (By default: \code{accession = NULL})
#'
#' @param as.df A single boolean value specifying whether to return the
#' results as a data.frame or leave as a nested list. In most cases,
#' \code{as.df = TRUE} will make the most sense.
#' (By default: \code{as.df = TRUE})
#'
#' @param max.hits A single integer value specifying the maximum number of
#' results to return or FALSE. The actual number of results will actually be
#' higher than \code{max.hits}, as clipping only occurs on pagination page
#' boundaries. To disable the limit, set \code{max.hits = NULL}.
#' (By default: \code{max.hits = 200})
#'
#' @param use.cache A single boolean value specifying whether to cache the
#' result -- and reuse any existing cache entry instead of issuing a new
#' callout. In general the use of caching for queries is discouraged, as
#' new data is being uploaded to MGnify all the time, which might potentially
#' be missed. However, for some purposes (such as analysis reproducibility)
#' caching makes sense. (By default: \code{use.cache = FALSE})
#'
#' @param ... Remaining parameter key/value pairs may be supplied to filter
#' the returned values. Available options differ between \code{types}.
#' See discussion above for details.
#'
#' @return A nested list or data.frame containing the results of the query.
#'
#' @examples
#' \dontrun{
#' mg <- MgnifyClient(cache_dir="/tmp/mgcache")
#'
#' # Get a list of studies from the Agricultural Wastewater :
#' agwaste_studies <- doQuery(
#'     mg, "studies", biome_name="Agricultural wastewater"
#'     )
#'
#' # Get all samples from a particular study
#' samps <- doQuery(mg, "samples", study_accession="MGYS00004521")
#'
#' # Search for all polar samples
#' samps_np <- doQuery(mg, "samples", latitude_gte=66, max.hits=NULL)
#' samps_sp <- doQuery(mg, "samples", latitude_lte=-66, max.hits=NULL)
#' samps_polar <- rbind(samps_np, samps_sp)
#' }
#'
#' @name doQuery
NULL

#' @rdname doQuery
#' @include MgnifyClient.R
#' @importFrom dplyr bind_rows
#' @export
setGeneric("doQuery", signature = c("x"), function(
        x, type = c("studies", "samples", "runs", "analyses"),
        accession = NULL, as.df = TRUE, max.hits = 200, use.cache = FALSE,
        ...
        )
    standardGeneric("doQuery"))

#' @rdname doQuery
#' @export
setMethod("doQuery", signature = c(x = "MgnifyClient"), function(
        x, type = c("studies", "samples", "runs", "analyses"),
        accession = NULL, as.df = TRUE, max.hits = 200, use.cache = FALSE,
        ...
        ){
    ############################### INPUT CHECK ################################
    if( !(.is_non_empty_string(type)) ){
        stop("'type' must be a single character value specifying ",
             "the type of instance to query.", call. = FALSE)
    }
    type <- match.arg(type, several.ok = FALSE)
    if( !(.is_non_empty_character(accession) || is.null(accession)) ){
        stop("'accession' must be a single character value or list of ",
             "character values specifying the MGnify accession identifier ",
             "or NULL.",
             call. = FALSE)
    }
    if( !.is_a_bool(as.df) ){
        stop("'as.df' must be a single boolean value specifying whether",
             "to return list or data.frame.", call. = FALSE)
    }
    if( !((.is_an_integer(max.hits) && (max.hits > 0 || max.hits == -1) ) ||
          is.null(max.hits) )  ){
        stop("'max.hits' must be a single integer value specifying the maximum ",
             "number of results to return or NULL.", call. = FALSE)
    }
    if( !.is_a_bool(use.cache) ){
        stop("'use.cache' must be a single boolean value specifying whether to ",
             "use on-disk caching.", call. = FALSE)
    }
    ############################# INPUT CHECK END ##############################
    # Perform query
    result <- .perform_query(client = x, type = type, accession = accession,
                             max.hits = max.hits, use.cache = use.cache, ...)
    # Convert list to data.frame if specified
    if(as.df){
        result <- .list_to_dataframe(result)
    }
    return(result)
})

################################ HELP FUNCTIONS ################################

.perform_query <- function(client, type, accession, max.hits, use.cache, ...){
    # Get optional arguments that were passed with ...
    qopt_list <- c(list(...), accession=accession)
    # Combine all arguments together
    all_query_params <- unlist(list(c(list(
        client = client,
        max.hits = max.hits,
        path = type,
        use.cache = use.cache,
        qopts = qopt_list
    ))), recursive = FALSE)
    # Get results from the database
    result <- do.call(".mgnify_retrieve_json", all_query_params)

    # Rename entries by accession
    id_list <- lapply(result, function(res) res$id)
    if( !is.null(result) ){
        names(result) <- id_list
    }
    return(result)
}

.list_to_dataframe <- function(result){
    dflist <- list()
    # Because metadata might not match across studies, the full dataframe
    # is built by first building per-sample dataframes, then using
    # rbind. fill from plyr to combine. For most use cases the number of
    # empty columns will hopefully be minimal... because who's going to
    # want cross study grabbing (?)
    for(r in result){
        df2 <- .mgnify_attr_list_to_df_row(
            json = r, metadata_key = "sample-metadata")
        # Loop through different datasets (e.g., biomes) that are related
        # to data
        for(rn in 1:length(r$relationships)){
            # Get specific relationship
            temp = r$relationships[[rn]]
            # Get only data of it
            temp_data <- temp$data
            # If there is data, include it
            # names(temp_data) %in% "id"
            if( !is.null(temp_data) ){
                # Take all "id" values. Some data can also include list of
                # lists. --> unlist and take "id" values
                temp_data <- unlist(temp_data)
                temp_data <- temp_data[names(temp_data) %in% "id"]
                temp_names <- rep(names(r$relationships)[rn],
                                  length(temp_data))
                # Get all column names and make them unique
                colnames <- append(colnames(df2), temp_names)
                colnames <- make.unique(colnames)
                # Get only column values that are being added
                temp_names <- colnames[
                    (length(colnames)-length(temp_names)+1):
                        length(colnames)]
                # Add new data to dataset
                df2[temp_names] <- temp_data
            }
        }
        # Add type of data that is being queried and accession code
        df2$type <- r$type
        rownames(df2) <- df2$accession
        # Add data to list
        dflist[[df2$accession]] <- df2
    }
    # Combine all data frames together
    result = bind_rows(dflist)
    return(result)
}
# # REMOVE THESE SINCE THEY ARE NOT USED?
# # Combined together into a single queriably list
# query_filters <- list(
#     biomes <- biome_filters,
#     samples <- sample_filters,
#     studies <- study_filters,
#     runs <- run_filters
# )
#
# # Filters possible - this comes from the django source code - would be nice if we could
# # look it up.
# # These DON'T seem to include all possible attributes ....
# # And only some
# sample_filters <- c(
#     'accession', 'experiment_type', 'biome_name', 'lineage', 'geo_loc_name',
#     'latitude_gte', 'latitude_lte', 'longitude_gte', 'longitude_lte',
#     'species','instrument_model', 'instrument_platform', 'metadata_key',
#     'metadata_value_gte', 'metadata_value_lte', 'metadata_value',
#     'environment_material', 'environment_feature', 'study_accession', 'include')
# biome_filters <- c(
#     'depth_gte', 'depth_lte')
# study_filters <- c(
#     'accession', 'biome_name', 'lineage', 'centre_name', 'include')
# run_filters <- c(
#     'accession', 'experiment_type', 'biome_name', 'lineage' , 'species',
#     'instrument_platform','instrument_model',
#     # 'metadata_key', 'metadata_value_gte', 'metadata_value_lte',
#     # 'metadata_value','sample_accession','study_accession',
#     'include')
# analysis_filters <- c(
#     'biome_name', 'lineage', 'experiment_type', 'species', 'sample_accession',
#     'pipeline_version')
