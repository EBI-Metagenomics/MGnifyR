#' Search MGnify database for studies, samples, runs, analyses, biomes,
#' assemblies, and genomes.
#'
#' @details
#' \code{doQuery} is a flexible query function, harnessing the "full"
#' power of the JSONAPI MGnify search filters. Search results may be filtered
#' by metadata value, associated study/sample/analyse etc.
#'
#' See [Api browser](https://www.ebi.ac.uk/metagenomics/api/v1/) for
#' information on MGnify database filters.
#' You can find help on customizing queries from
#' [here](https://emg-docs.readthedocs.io/en/latest/api.html#customising-queries).
#'
#' For example the following filters are available:
#' \itemize{
#'     \item{\strong{studies}: accession, biome_name, lineage, centre_name,
#'     include}
#'     \item{\strong{samples}: accession, experiment_type, biome_name,
#'     lineage, geo_loc_name, latitude_gte, latitude_lte,
#'     longitude_gte, longitude_lte, species, instrument_model,
#'     instrument_platform, metadata_key, metadata_value_gte,
#'     metadata_value_lte, metadata_value, environment_material,
#'     environment_feature, study_accession, include}
#'     \item{\strong{runs}: accession, experiment_type, biome_name, lineage,
#'     species, instrument_platform, instrument_model, metdata_key,
#'     metadata_value_gte, metadata_value_lte, metadata_value, sample_accession,
#'     study_accession, include}
#'     \item{\strong{analyses}: biome_name, lineage, experiment_type, species,
#'     sample_accession, pipeline_version}
#'     \item{\strong{biomes}: depth_gte, depth_lte}
#'     \item{\strong{assemblies}: depth_gte, depth_lte}
#'  }
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
#' \code{runs}, \code{analyses}, \code{biomes}, \code{assemblies},
#' \code{super-studies}, \code{experiment-types}, \code{pipelines},
#' \code{pipeline-tools}, \code{publications}, \code{genomes},
#' \code{genome-search}, \code{genome-search/gather}, \code{genome-catalogues},
#' \code{genomeset}, \code{cogs}, \code{kegg-modules}, \code{kegg-classes},
#' \code{antismash-geneclusters}, \code{annotations/go-terms},
#' \code{annotations/interpro-identifiers}, \code{annotations/kegg-modules},
#' \code{annotations/pfam-entries}, \code{annotations/kegg-orthologs},
#' \code{annotations/genome-properties},
#' \code{annotations/antismash-gene-clusters}, \code{annotations/organisms}, or
#' \code{mydata}.
#' (By default: \code{type = "studies"})
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
#' @param ... Remaining parameter key/value pairs may be supplied to filter
#' the returned values. Available options differ between \code{types}.
#' See discussion above for details.
#'
#' @return A nested list or data.frame containing the results of the query.
#'
#' @examples
#' mg <- MgnifyClient(useCache = FALSE)
#'
#' # Get a list of studies from the Agricultural Wastewater :
#' agwaste_studies <- doQuery(
#'     mg, "studies", biome_name="Agricultural wastewater"
#'     )
#'
#' \dontrun{
#' # Get all samples from a particular study
#' samps <- doQuery(mg, "samples", accession="MGYS00004521")
#'
#' # Search polar samples
#' samps_np <- doQuery(mg, "samples", latitude_gte=66, max.hits=10)
#' samps_sp <- doQuery(mg, "samples", latitude_lte=-66, max.hits=10)
#'
#' # Search studies that have studied drinking water
#' tbl <- doQuery(
#'     mg,
#'     type = "studies",
#'     biome_name = "root:Environmental:Aquatic:Freshwater:Drinking water",
#'     max.hits = 10)
#' }
#'
#' @name doQuery
NULL

#' @rdname doQuery
#' @importFrom dplyr bind_rows
#' @include AllClasses.R AllGenerics.R MgnifyClient.R utils.R
#' @export
setMethod("doQuery", signature = c(x = "MgnifyClient"), function(
        x, type = "studies", accession = NULL, as.df = TRUE, max.hits = 200,
        ...){
    ############################### INPUT CHECK ################################
    available_types <- c(
        "studies", "samples", "runs", "analyses", "biomes", "assemblies",
        "super-studies", "experiment-types", "pipelines", "pipeline-tools",
        "publications", "genomes", "genome-search", "genome-search/gather",
        "genome-catalogues", "genomeset", "cogs", "kegg-modules",
        "kegg-classes", "antismash-geneclusters", "annotations/go-terms",
        "annotations/interpro-identifiers", "annotations/kegg-modules",
        "annotations/pfam-entries", "annotations/kegg-orthologs",
        "annotations/genome-properties", "annotations/antismash-gene-clusters",
        "annotations/organisms", "mydata")
    if( !(.is_non_empty_string(type) && type %in% available_types) ){
        stop(
            "'type' must be a single character value specifying ",
            "the type of instance to query. The value must be one of the ",
            "following options: ",
            paste0("'", paste(available_types, collapse = "', '"), "'"),
            call. = FALSE)
    }
    if( !(.is_non_empty_character(accession) || is.null(accession)) ){
        stop(
            "'accession' must be a single character value or vector of ",
            "character values specifying the MGnify accession identifier ",
            "or NULL.",
            call. = FALSE)
    }
    if( !.is_a_bool(as.df) ){
        stop(
            "'as.df' must be a single boolean value specifying whether",
            "to return list or data.frame.", call. = FALSE)
    }
    if( !((.is_an_integer(max.hits) && (max.hits > 0 || max.hits == -1) ) ||
        is.null(max.hits) )  ){
        stop(
            "'max.hits' must be a single integer value specifying the ",
            "maximum number of results to return or NULL.", call. = FALSE)
    }
    ############################# INPUT CHECK END ##############################
    # Perform query
    result <- .perform_query(
        client = x, type = type, accession = accession, max.hits = max.hits,
        ...)
    # Convert list to data.frame if specified
    if( as.df && length(result) > 0 ){
        # Turn lists to dfs
        result <- lapply(result, .list_to_dataframe)
        # Combine dfs
        result <- bind_rows(result)
    }
    return(result)
})

################################ HELP FUNCTIONS ################################

.perform_query <- function(
        client, type, accession, max.hits, use.cache = useCache(client),
        show.messages = verbose(client), ...){
    # Input check
    if( !.is_a_bool(use.cache) ){
        stop(
            "'use.cache' must be a single boolean value.", call. = FALSE)
    }
    #
    # Get parameters that are passed to do the query from database
    query_params <- list(...)
    query_params[["accession"]] <- accession
    # Get results from the database
    result <- .mgnify_retrieve_json(
        client = client,
        path = type,
        max.hits = max.hits,
        use.cache = use.cache,
        qopts = query_params
        )
    # Rename entries by accession
    id_list <- lapply(result, function(res) res$id)
    if( !is.null(result) ){
        names(result) <- id_list
    }
    return(result)
}

.list_to_dataframe <- function(result){
    # Get attributes
    df <- .mgnify_attr_list_to_df_row(
      json = result, metadata_key = "sample-metadata")

    # Loop through relationships, i.e., this data might be related to specific
    # samples, analyses... --> get that info
    relationships <- result[["relationships"]]
    for( i in seq_len(length(relationships)) ){
        # Get specific relationship, e.g., this data vs related runs
        relationship_type <- names(result$relationships)[[i]]
        relationship <- result$relationships[[i]]
        # Get only data (temp is list of lists and only data element contains
        # relevant info)
        rel_data <- relationship[["data"]]
        # If there is data, include it
        if( !is.null(rel_data) && length(rel_data) > 0 ){
            # Take all "id" values. Some data can also include list of
            # lists. --> unlist and take "id" values. Based on this ID (such
            # as "runs" ID) user can fetch specific relationship
            rel_data <- unlist(rel_data)
            rel_data <- rel_data[names(rel_data) %in% "id"]
            temp_names <- rep(relationship_type, length(rel_data))
            # Get all column names and make them unique
            colnames <- append(colnames(df), temp_names)
            colnames <- make.unique(colnames)
            # Get only column values that are being added
            temp_names <- colnames[
                (length(colnames)-length(temp_names)+1):length(colnames)]
            # Add new data to dataset
            df[temp_names] <- rel_data
        }
    }
    # Add type of data that is being queried and accession code
    df[["type"]] <- result[["type"]]
    rownames(df) <- df[["accession"]]
    return(df)
}
