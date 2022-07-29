#' Search MGnify database for studies, samples and runs
#'
#' \code{mgnify_query} is a flexible query function, harnessing the "full" power of the JSONAPI MGnify
#' search filters. Search results may be filtered by metadata value, associated study/sample/analyese etc. Details
#' of the capabilities may be found [here](https://emg-docs.readthedocs.io/en/latest/api.html#customising-queries). Currently,
#' the following filters are available (based on examination of the Python source code):
#'     \itemize{
#'         \item{\strong{Studies} : accession, biome_name, lineage, centre_name}
#'            \item{\strong{Samples} : accession, experiment_type, biome_name,
#'     lineage, geo_loc_name, latitude_gte, latitude_lte,
#'     longitude_gte, longitude_lte, species, instrument_model, instrument_platform,
#'        metadata_key, metadata_value_gte, metadata_value_lte, metadata_value,
#'        environment_material, environment_feature, study_accession}
#'        \item{\strong{Runs} accession, experiment_type, biome_name, lineage, species,
#'            instrument_platform, instrument_model}
#'        }
#'        Unfortunately it appears that in some cases, some of these filters don't work as expected, so it is important
#'        to check the results returned match up with what's expected. Even more unfortunately if there's an error in the
#'        parameter specification, the query will run as if no filter parameters were present at all. Thus the
#'        result will appear superficially correct but will infact correspond to something completely different. This behaviour
#'        will hopefully be fixed in future incarnations of the MGnifyR or JSONAPI, but for now users should double check returned
#'        values.
#'
#'        It is currently not possible to combine queries of the same type in a single call (for example to search for samples
#'        \emph{between} latitude). However, it is possible to run multiple queries and combine the results using set operations in R to get the
#'        desired behaviour.
#'
#'
#' @importFrom dplyr bind_rows
#'
#' @param client (To be described)
#' @param qtype Type of objects to query. One of \code{studies},\code{samples},\code{runs} or
#' \code{analyses}
#' @param accession Either a single known MGnify accession identifier (of type \code{qtype}), or a list/vector
#' of accessions to query. Note that multiple values only work for samples, runs and assemblies ... not sure why.
#' @param asDataFrame Boolean flag to choose whether to return the results as a data.frame or leave as a nested list. In
#' most cases, \code{asDataFrame = TRUE} will make the most sense.
#' @param maxhits determines the maximum number of results to return. The actual number of results will actually be higher than \code{maxhits},
#' as clipping only occurs on pagination page boundaries. To disable the limit, set \code{maxhits} < 0
#' @param usecache Whether to cache the result - and reuse any existing cache entry instead of issuing a
#' new callout. In generl the use of caching for queries is discouraged, as new data is being uploaded to MGnify
#' all the time, which might potentially be missed. However, for some purposes (such as analysis reproducibility)
#' caching makes sense.
#' @param ... Remaining parameter key/value pairs may be supplied to filter the returned values. Available options differ
#' between \code{qtypes}.See discussion above for details.
#' @return A nested list or data.frame (depending on \code{asDataFrame}) containing the results of the query.
#' @examples
#' mg <- mgnify_client(cache_dir="/tmp/mgcache")
#'
#' #Get a list of studies from the Agricultural Wastewater :
#' agwaste_studies <- mgnify_query(mg, "studies", biome_name="Agricultural wastewater")
#'
#' #Get all samples from a particular study
#' samps <- mgnify_query(mg, "samples", study_accession="MGYS00004521")
#'
#' #Search for all polar samples
#' samps_np <- mgnify_query(mg, "samples", latitude_gte=66, maxhits=-1)
#' samps_sp <- mgnify_query(mg, "samples", latitude_lte=-66, maxhits=-1)
#' samps_polar <- rbind(samps_np, samps_sp)
#'
#' @export
mgnify_query <- function(client, qtype="samples", accession=NULL, asDataFrame=T, maxhits=200, usecache=F, ...){
    #Need to get around the lazy expansion in R in order to get a list
    a <- accession
    arglist <- as.list(match.call())[-1] # drop off the first entry, which is the name of the function

    arglist$accession <- a
    #Force evaluation of arguments to prevent things getting messed up with x and y and z and ....
    ## arglist <- lapply(arglist, eval)
    #Filter the query options such that
    #qopt_list <- arglist[names(arglist) %in% query_filters[[qtype]]]
    qopt_list <- c(list(...), accession=accession)
    #qopt_list <- lapply(qopt_list, force)
    non_qopts <- arglist[!(names(arglist) %in% c(c("asDataFrame","qtype","client", "maxhits"),query_filters[[qtype]]))]

    all_query_params <- unlist(list(c(list(client=client, maxhits=maxhits, path=qtype, usecache=usecache, qopts=qopt_list))), recursive = F)

    result <- do.call("mgnify_retrieve_json", all_query_params)

    #Rename entries by accession
    id_list <- lapply(result, function(x) x$id)
    names(result) <- id_list

    if(asDataFrame){
        #Because metadata might not match across studies, the full dataframe is built by first building per-sample dataframes,
        #then using rbind.fill from plyr to combine. For ~most~ use cases the number of empty columns will hopefully
        #be minimal... because who's going to want cross study grabbing (?)
        dflist <- lapply(result, function(r){
            df2 <- mgnify_attr_list_to_df_row(json = r, metadata_key = "sample-metadata")

            #Currently this is a bit hacky - assumes the study only has one biome, and sample only one study etc.
            for(rn in names(r$relationships)){
                tryCatch({
                    #I truely appologise for the following line of code. It's AWFUL! But it handles
                    #the case where data is both a list or a list of lists... I think.
                    df2[rn] <- unlist(list(a=list(r$relationships[[rn]]$data),b=list()), recursive = F, use.names = F)[[1]]$id
                },
                error <- function(x)warning(paste("Failed to add entry for ", rn, " to ", df2$accession[[1]], sep="")))
            }
            df2$type <- r$type
            rownames(df2) <- df2$accession
            df2
        }
        )
        tryCatch(
            dplyr::bind_rows(dflist),
            error <- function(e) dflist
        )
    }else{
        result
    }

}
