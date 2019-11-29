# NOTES:

#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


library(httr)
library(urltools)
library(phyloseq)
#library(plyr)
library(reshape2)
library(dplyr)



##Example accessions:
# Analysis assembly: MGYA00379728
# Analysis metagenomic: MGYA00377678
# Analysis amplicon: MGYA00250889



baseurl='https://www.ebi.ac.uk/metagenomics/api/v1'

##Filters possible - this comes from the django source code - would be nice if we could
# look it up.
# These DON'T seem to include all possible attributes ....
# And only some
sample_filters = c('accession','experiment_type','biome_name','lineage','geo_loc_name','latitude_gte','latitude_lte',
                   'longitude_gte','longitude_lte','species','instrument_model','instrument_platform','metadata_key',
                   'metadata_value_gte','metadata_value_lte','metadata_value','environment_material','environment_feature',
                   'study_accession','include')
biome_filters = c('depth_gte','depth_lte')
study_filters = c('accession','biome_name','lineage','centre_name','include')
run_filters=c('accession','experiment_type','biome_name','lineage','species','instrument_platform','instrument_model',
              # 'metadata_key','metadata_value_gte','metadata_value_lte','metadata_value','sample_accession','study_accession',
              'include')
analysis_filters = c('biome_name', 'lineage', 'experiment_type', 'species', 'sample_accession', 'pipeline_version')


#Combined together into a single queriably list
query_filters=list(
  biomes=biome_filters,
  samples=sample_filters,
  studies  = study_filters,
  runs=run_filters
)





#**************************************************
#Classes

#' MGnify API client object.
#'
#' Acts as a simple container encapsulating API info (\code{baseurl}), user parameters (\code{authtok}) and cache
#' locations (\code{cache_dir})
#'
#' @slot baseurl Web address of the MGnify JSON api - including version. Required
#' @slot authtok The MGnify API supports authentication for users by way of \code{authtoken}. This is that token, althouth it's not used in anger at the moment.
#' @slot cache_dir To reduce load on the server, and speed up repeated data processing, JSON calls may be cached locally, and
#' stored in \code{cache_dir}.
#' @export mgnify_client
##' @exportClass mgnify_client
mgnify_client <- setClass("mgnify_client",
         slots=list(url = "character", authtok = "character", cache_dir="character"),
         prototype = list(url=baseurl, authtok=NULL, cache_dir=NULL))

#Contructor to allow logging in with username/password
mgnify_client <- function(username=NULL,password=NULL,usecache=F,cachedir=NULL){
  url=baseurl
  authtok=NA_character_

  #Check to see if we're goint to try and get an authentication token:
  if (!is.null(username) & !is.null(password)){
    r = httr::POST(paste(url, "utils/token/obtain", sep="/"),
                   body=list(username=username, password=password),
                   encode="json")
    cont = httr::content(r)
    if ("data" %in% names(cont)){
      authtok = cont$data$token
    }
    else{
      "Failed to authenticate"
    }
  }
  #Assume we're not using it
  cachepath=NULL
  if(usecache){
    if (is.null(cachedir) ){
      cachepath=paste(getwd(),'.MGnifyR_cache',sep="/")
    }else{
      cachepath=cachedir
    }
    #Make it if needed - assume the user is sensible and the path will work...
    dir.create(cachepath,showWarnings = F)
  }

  #Return the final object
  new("mgnify_client", url=url, authtok=authtok, cache_dir = cachepath)

}



#helper function for getting relative paths in the API
#Not everything is implemented here - just what we
#need to get to the download or run areas
#Given an accession x, we want to get the link to get the url for the
#coresponding typeY
#' \code{JSONAPI} path for child elements
#'
#' \code{mgnify_get_x_for_y} determines the location of \code{typeY} child objects of \code{x} (type \code{typeX})
#'
#' This helper function, principally intended to be used internally, is used to match up related objects within the path. The inherently
#' unhierarchical nature of the MGnify API makes it a bit inconsistent. This function acts as a quick way to determine how to get from
#' one type to another, without having to special case within the code.
#' @param client MGnifyR client API object
#' @param x Accession ID \code{char} of parent object
#' @param typeX Type of accession \code{x}
#' @param typeY Type of child object to return
#' @return \code{char} complete url to access the result. Note this query is not run from here - just the URL is returned
#' @examples
#' cl <- new("mgnify_client")
#' mgnify_get_x_for_y(cl, "MGYS00005126", "studies", "samples")
##'@export
mgnify_get_x_for_y <- function(client, x, typeX, typeY, usecache=F){
  #This one's easy - just rearrange the URLs
  #if(typeX=="samples" & typeY %in% c("runs","studies")){
  #  paste( typeX,x,typeY, sep="/")
  #}else if(typeX=="runs" & typeY == "analyses"){
  #  paste( typeX,x,typeY, sep="/")
  #}
  #else{
    #Do it the hard way with a callout
    json_dat = mgnify_retrieve_json(client, paste(typeX, x, sep="/"), usecache = usecache)
    #cat(str(json_dat))
    #tgt_access = json_dat[[1]]$relationships[[typeY]]$data$id
    #tgt_type = json_dat[[1]]$relationships[[typeY]]$data$type
    #paste(tgt_type,tgt_access,sep="/")
    json_dat[[1]]$relationships[[typeY]]$links$related
    #substr(tgt_url, nchar(client@url) + 1, nchar(tgt_url))
  #}
}


#' Look up analysis accession IDs for one or more study accessions
#'
#' \code{mgnify_analyses_from_studies} Retrieve Analysis accession IDs associated with the supplied Study accession
#'
#' Helper function to get all analyses associated with the given studies.
#'
#' @param client \code{mgnify_client} instance
#' @param accession Single study accession id, or vector/list of accessions for which to retrieve Analyses ids
#' @param usecache Flag to determine whether to re-use/store data on disk, rather than query the server.
#' @return vector of Analysis accession ids
#'
#' @examples
#' #Retrieve all analysis ids from studies MGYS00005058, MGYS00005058 and MGYS00005058
#' result <- mgnify_analyses_from_studies(myclient, c("MGYS00005058", "MGYS00005058" and "MGYS00005058"))
#'
#' @export
mgnify_analyses_from_studies <- function(client, accession, usecache=F){
    analyses_accessions <- sapply(as.list(accession), function(x){
      accurl <- mgnify_get_x_for_y(client, x, "studies","analyses", usecache = usecache )
      jsondat <- mgnify_retrieve_json(client, complete_url = accurl, usecache = usecache, maxhits = -1)
      #Just need the accession ID
      lapply(jsondat, function(x) x$id)
    })
    unlist(analyses_accessions)
}


#' Look up analysis accession IDs for one or more sample accessions
#'
#' \code{mgnify_analyses_from_samples} Retrieve Analysis accession IDs associated with the supplied Sample accession
#'
#' Helper function to get all analyses associated with the given samples.
#'
#' @param client \code{mgnify_client} instance
#' @param accession Single sample accession id, or vector/list of accessions for which to retrieve Analyses ids
#' @param usecache Flag to determine whether to re-use/store data on disk, rather than query the server.
#' @return vector of associated Analysis accession ids
#'
#' @examples
#' #Retrieve all analysis ids from samples
#' result <- mgnify_analyses_from_studies(myclient, c("MGYS00005058", "MGYS00005058" and "MGYS00005058"))
#'
#' @export
mgnify_analyses_from_samples <- function(client, accession, usecache=F){
  analyses_accessions <- sapply(as.list(accession), function(x){
    accurl <- mgnify_get_x_for_y(client, x, "samples","analyses", usecache = usecache )
    jsondat <- mgnify_retrieve_json(client, complete_url = accurl, usecache = usecache)
    #Just need the accession ID
    lapply(jsondat, function(x) x$id)
  })
  unlist(analyses_accessions)
}



#Not exporting this - if people want to they can use the
# rjsonapi functionality
#'Coverting attribute lists to a single data.frame row
#'
#'\code{mgnify_attr_list_to_df_row} extracts the \code{attribute} entry in a JSONAPI data result
#'and converts it to a single row data.frame with columns corresponding to attribute entries
#'Optionally parses a \code{metadata_key} item for extra metadata items.
#'@param json The \emph{raw} result list
#'@param metadata_key Optional extra key to parse subattributes from
#'@return data.frame containing a single row of metadata
#'
#'@export
mgnify_attr_list_to_df_row <- function (json, metadata_key=NULL ){

  attrlist=names(json$attributes)

  if (!is.null(metadata_key)){
    baseattrlist=attrlist[!(attrlist %in% c(metadata_key))]
    metaattrlist=json$attributes[[metadata_key]]
    metlist=sapply(metaattrlist, function(x) x$value)
    names(metlist)=sapply(metaattrlist, function(x) x$key)
    df = as.data.frame(t(unlist(c(json$attributes[baseattrlist], metlist))), stringsAsFactors = F)
  }else{
    df = as.data.frame(t(unlist(json["attributes"])), stringsAsFactors = F)
  }
  df$accession <- json$id
  df$acc_type <- json$type

  rownames(df)=df$accession
  df
}



#Internal function to actually perform the http request. Build up the URL then issues
#a GET, parsing the returned JSON into a nested list (uses \code{jsonlite} internally?)
#Previously cached results may be retrieved from disk without resorting to calling the MGnify server.

#'Low level MGnify API handler
#'
#'\code{mgnify_retrieve_json} deals with handles the actual HTTP GET calls for the MGnifyR package, handling both pagination and local reuslt
#'caching. Although principally intended for internal MGnifyR use , it's exported for direct invocation.
#'
#'@param client MGnifyR client
#'@param path top level search point for the query. One of \code{biomes}, \code{samples}, \code{runs} etc.
#'@param complete_url \emph{complete} url to search, usuaally retrieved from previous query in the "related" section.
#'@param qopts named list or vector containing options/filters to be URL encoded and appended to query as key/value pairs
#'@param maxhits Maxmium number of data entries to return. The actual number of hits returned may be higher than this value,
#'as this parameter only clamps after each full page is processed. Set to <=0 to disable - i.e. retrieve all items.
#'@param usecache Should successful queries be cached on disk locally? There are unresolved questions about whether this is
#'a sensible thing to do, but it remains as an option. It probably makes sense for single accession grabs, but not for
#'(filtered) queries - which are liable to change as new data is added to MGnify. Also caching only works for the first page.
#'@param Debug Should we print out lots of information while doing the grabbing?
#'@return \code{list} of results after pagination is dealt with.
#'@export
  mgnify_retrieve_json <- function(client, path="biomes", complete_url=NULL, qopts=NULL,
                                   maxhits=200, usecache = F, Debug=F){
  #Set up the base url
  #Are we using internal paths?
  if (is.null(complete_url)){
    fullurl = paste(client@url, path, sep="/")
  }
  #Or direct links from e.g. a "related" section
  else{
    #Set the full url, but clean off any existing parameters (page, format etc) as they'll be added back later:
    fullurl = complete_url
    urltools::parameters(fullurl) <- NULL
    path = substr(fullurl, nchar(client@url) + 2, nchar(fullurl))
  }

  #cat(fullurl)

  #convert to csv if filters are lists.
  #This doesn't check if they ~can~ be searched for in the API,
  #which is an issue since no error is returned by the JSON if the search
  #is invalid - we only get a result as if no query was present...
  tmpqopts = lapply(qopts,function(x) paste(x,collapse = ','))

  #Include the json and page position options
  #full_qopts = as.list(c(format="json", tmpqopts, page=1))
  full_qopts = as.list(c(format="json", tmpqopts))
  #Build up the cache name anyway - even if it's not ultimately used:
  fname_list = c(path, names(unlist(full_qopts)), unlist(full_qopts))
  cache_fname = paste(fname_list,collapse = "_")
  cache_full_fname = paste(client@cache_dir, '/', cache_fname, '.RDS', sep="")
  # Do we want to try and use a cache to speed things up?

  if(usecache & file.exists(cache_full_fname)){
      final_data = readRDS(cache_full_fname)
  }else{

    res = httr::GET(url=fullurl, config(verbose=Debug), query=full_qopts )
    data <-httr::content(res)

    #At this point, data$data is either a list of lists or a single named list. If it's a single entry, it needs embedding in
    #a list for consistency downstream
    #datlist is built up as a list of pages, where each entry must be another list. Thus, on the first page,
    #
    datlist=list()
    if (!is.null(names(data$data))){
    #Create something to store the returned data

      datlist[[1]] = list(data$data)
    }else{
      datlist[[1]] = data$data
    }
      #cat(str(data))
    # Check to see if there's pagination required
    if ("meta" %in% names(data)){
      #Yes, paginate
      pstart = as.numeric(data$meta$pagination$page)
      pend   = as.numeric(data$meta$pagination$pages)

      for (p in seq(pstart+1,pend)){  # We've already got the first one

        full_qopts$page=p
        curd = httr::content(httr::GET(fullurl, config(verbose=Debug), query=full_qopts ))
        datlist[[p]] = curd$data
        #Check to see if we've pulled enough entries
        if(maxhits > 0){
          curlen=sum(sapply(datlist, length))
          if (curlen > maxhits){
            break
          }
        }
      }
    }
    #if(length(datlist) > 1){
    final_data <- unlist(datlist, recursive=F)
    #}else{
    #  final_data <- datlist
    #}
    #}
    #else{
    #  final_data <- datlist
    if (usecache && !file.exists(cache_full_fname)){
      #Make sure the directory is created...
      dir.create(dirname(cache_full_fname), recursive = T)
      saveRDS(final_data, file = cache_full_fname)
    }
  }
  final_data
}


#'#' Search MGnify database for studies, samples and runs
#'
#' \code{mgnify_query} is a flexible query function, harnessing the "full" power of the JSONAPI MGnify
#' search filters.
#' @param mgnify_client instance
#' @param qtype Type of objects to query. One of \code{studies},\code{samples},\code{runs} or
#' \code{analyses}
#' @param accession Either a single known MGnify accession identifier (of type \code{qtype}), or a list/vector
#' of accessions to query.
#' @param asDataFrame Boolean flag to choose whether to return the results as a data.frame or leave as a nested list. In
#' most cases, \code{asDataFrame = TRUE} will make the most sense.
#' @param maxhits determines the maximum number of results to return. The actual number of results will actually be higher than \code{maxhits},
#' as clipping only occurs on pagination page boundaries. To disable the limit, set \code{maxhits} < 0
#' @export
mgnify_query <- function(client, qtype="samples", accession=NULL, asDataFrame=F, maxhits=200, ...){
  #Need to get around the lazy expansion in R in order to get a list
  a=accession
  arglist = as.list(match.call())[-1] # drop off the first entry, which is the name of the function
#  arglist

  arglist$accession=a

  #Filter the query options such that
  qopt_list = arglist[names(arglist) %in% query_filters[[qtype]]]
  non_qopts = arglist[!(names(arglist) %in% c(c("asDataFrame","qtype","client", "maxhits"),query_filters[[qtype]]))]

  #cat(str(arglist))
  all_query_params = unlist(list(c(list(client=client, maxhits=maxhits, path=qtype, qopts=qopt_list))), recursive = F)
  cat(str(all_query_params))
  #cat(str(all_query_params))
  #Do the query
  #result = mgnify_retrieve_json(client, path=qtype, qopts = qopts)
  result = do.call("mgnify_retrieve_json", all_query_params)

  #Rename entries by accession
  id_list = lapply(result, function(x) x$id)
  names(result) = id_list

  if(asDataFrame){
    #Because metadata might not match across studies, the full dataframe is built by first building per-sample dataframes,
    # then using rbind.fill from plyr to combine. For ~most~ use cases the number of empty columns will hopefully
    # be minimal... because who's going to want cross study grabbing (?)
    dflist = lapply(result, function(r){
      df2 <- mgnify_attr_list_to_df_row(json = r, metadata_key = "sample-metadata")
      df2$biome = r$relationship$biome$data$id
      df2$study = r$relationship$studies$data$id
      df2$type = r$type
      rownames(df2)=df2$accession
      df2

    }
    )
    tryCatch(
      dplyr::bind_rows(dflist),
      error=function(e) dflist
    )
  }else{
    result
  }

}


#Retrieves combined study/sample/analysis metadata
mgnify_get_single_analysis_metadata <- function(client=NULL, accession, usecache=T){

  dat <- mgnify_retrieve_json(client, paste("analyses", accession, sep="/"), usecache = usecache)
  #There ~should~ be just a single result
  top_data <- dat[[1]]
  analysis_df <- mgnify_attr_list_to_df_row(top_data, metadata_key = "analysis-summary")

  #cat(str(analysis_metadata))
  #Build up the metadata dataframe from the analyses_metadata_headers vector:
  sample_met <- mgnify_retrieve_json(client, complete_url = top_data$relationships$sample$links$related, usecache = usecache)
  study_met <- mgnify_retrieve_json(client, complete_url = top_data$relationships$study$links$related, usecache = usecache)

  sample_df <- mgnify_attr_list_to_df_row(sample_met[[1]], metadata_key = "sample-metadata")
  study_df <- mgnify_attr_list_to_df_row(study_met[[1]])

  colnames(sample_df) <- paste("sample",colnames(sample_df), sep="_")
  colnames(study_df) <- paste("study",colnames(study_df), sep="_")
  colnames(analysis_df) <- paste("analysis",colnames(analysis_df), sep="_")

  rownames(sample_df) <- rownames(analysis_df)
  rownames(study_df) <- rownames(analysis_df)
  cbind(analysis_df, study_df, sample_df)
}



#Internal function to parse the attributes/hierarchy list into a data.frame
mgnify_parse_tax <- function(json){
  df <- as.data.frame(c(json$attributes["count"], unlist(json$attributes$hierarchy)), stringsAsFactors = F)
  df$index_id <- json$attributes$lineage
  df

}

mgnify_parse_func <- function(json){
  df <- as.data.frame(json$attributes, stringsAsFactors = F)
  df$index_id <- json$attributes$accession
  df
}

analyses_results_type_parsers <- list(taxonomy=mgnify_parse_tax,`taxonomy-itsonedb` = mgnify_parse_tax, `go-slim`=mgnify_parse_func,
                                      `taxonomy-itsunite`=mgnify_parse_tax, `taxonomy-ssu`=mgnify_parse_tax,
                                      `taxonomy-lsu`=mgnify_parse_tax,`antismash-gene-clusters`=mgnify_parse_func,
                                      `go-terms`=mgnify_parse_func, `interpro-identifiers`=mgnify_parse_func)


#Retrieves combined study/sample/analysis metadata
mgnify_get_single_analysis_results <- function(client=NULL, accession, retrievelist=c(), usecache=T){
  metadata_df <- mgnify_get_single_analysis_metadata(client, accession, usecache=usecache)
  #Now (re)load the analysis data:
  analysis_data <- mgnify_retrieve_json(client, paste("analyses",accession,sep="/"), usecache = usecache)
  #For now try and grab them all - just return the list - don't do any processing...
  all_results <- lapply(names(analyses_results_type_parsers), function(r) {
    if(r %in% retrievelist){
      tmp <- mgnify_retrieve_json(client, complete_url = analysis_data[[1]]$relationships[[r]]$links$related, usecache = usecache)
    #cat(str(tmp))
    #cat(str(tmp))
      tmp
    }

  })
  names(all_results) <- names(analyses_results_type_parsers)

  parsed_results = sapply(names(all_results), function(x){
    all_json <- all_results[[x]]
    res_df <- do.call(dplyr::bind_rows, lapply(all_json,analyses_results_type_parsers[[x]] ))
    rownames(res_df) <- res_df$index_id
    res_df
  })
  parsed_results
}





mgnify_get_single_analysis_phyloseq <- function(client=NULL, accession, usecache=T, downloadDIR=NULL, tax_SU="SSU"){
  metadata_df <- mgnify_get_single_analysis_metadata(client, accession, usecache=usecache)

  analysis_data <-  mgnify_retrieve_json(client, paste("analyses",accession,sep="/"), usecache = usecache)
  download_url <- analysis_data[[1]]$relationships$downloads$links$related
  analysis_downloads <- mgnify_retrieve_json(client, complete_url = download_url,usecache = usecache)

  #Depending on the pipeline version, there may be more than one OTU table available (LSU/SSU), so try and get the
  #one specified in tax_SU - otherwise spit out a warning and grab the generic (older pipelines)
  available_biom_files <- analysis_downloads[grepl('JSON Biom', sapply(analysis_downloads, function(x) x$attributes$`file-format`$name))]
  biom_position <- grepl(tax_SU, sapply(available_biom_files, function(x) x$attributes$`group-type`))
  if(sum(biom_position) == 0){
    warning("Unable to locate requested taxonomy type ",tax_SU,". This is likely due to the current analysis having been performed on an older version of the MGnify pipeline.
             The available BIOM file will be used instead.")
    biom_url <- available_biom_files[[1]]$links$self
  }else{
    biom_url <- available_biom_files[biom_position][[1]]$links$self
  }

  #Can specify a seperate dir for saving biom files, otherwise they end up in the client@cachdir folder, under "bioms"
  if (is.null(downloadDIR)){
    downloadDIR=paste(client@cache_dir,"biom_files",sep="/")
    dir.create(downloadDIR, recursive = T, showWarnings = F)
  }
  fname=tail(strsplit(biom_url, '/')[[1]], n=1)
  biom_path = paste(downloadDIR, fname, sep="/")
  if (! file.exists(biom_path)){#} | !use_downloads ){
    httr::GET(biom_url, write_disk(biom_path, overwrite = T))
  }
  #Load in the phlyloseq object
  psobj <- phyloseq::import_biom(biom_path)
  #Need to check if the taxonomy was parsed correctly - depending on the pipeline it may need a bit of help:
  if (ncol(tax_table(psobj)) == 1){
    psobj <- phyloseq::import_biom(biom_path, parseFunction = parse_taxonomy_qiime)
  }
  if(! "Kingdom" %in% names(tax_table(psobj))){
    psobj <- phyloseq::import_biom(biom_path, parseFunction = parse_taxonomy_greengenes)
  }

  #The biom files have a single column of "sa1". It's rewritten as sample_run_analysis accession, with
  # the original value stored in the sample_data (just in case it changes between pipelines)
  orig_samp_name <- sample_names(psobj)[[1]]
  newsampname <- rownames(metadata_df)[1]
  metadata_df[1,"orig_samp_name"] <- orig_samp_name
  sample_names(psobj) <- newsampname
  sample_data(psobj) <- metadata_df
  psobj
}



mgnify_get_analysis_metadata <- function(client, accessions, usecache=T){
  reslist <- lapply(as.list(accessions), function(x) mgnify_get_single_analysis_metadata(client, x, usecache = T))
  df <- do.call(dplyr::bind_rows,reslist)
  rownames(df) <- accessions
  df
}

mgnify_get_analyses_phyloseq <- function(client = NULL, accessions, usecache=T, returnLists=F, tax_SU = "SSU"){
  #Some biom files don't import - so need a try/catch
  ps_list <- plyr::llply(accessions, function(x) {tryCatch(
                    mgnify_get_single_analysis_phyloseq(client, x, usecache = usecache, tax_SU = tax_SU), error=function(x) NULL)

    }, .progress = "text")

  #The sample_data has been corrupted by doing the merge (names get messed up and duplicated), so just regrab it with another lapply/rbind
  samp_dat <- lapply(accessions, function(x) mgnify_get_single_analysis_metadata(client, x, usecache = usecache ))
  if (returnLists){
    list(phyloseq_objects=ps_list, sample_metadata = samp_dat)
  }else{

    full_ps <- do.call(merge_phyloseq, ps_list)
    sample_metadata_df <- do.call(dplyr::bind_rows, samp_dat)
    rownames(sample_metadata_df) <- sample_metadata_df$analysis_accession
    sample_data(full_ps) <- sample_metadata_df
    full_ps
  }
}



mgnify_get_analysis_results <- function(client=NULL, accessions, retrievelist=c(), compact_results=T, usecache = T){
  if(length(retrievelist) == 1 && retrievelist == "all"){
    retrievelist = names(analyses_results_type_parsers)
  }
  results_as_lists <- plyr::llply(accessions, function(x) mgnify_get_single_analysis_results(client, x, usecache = usecache, retrievelist = retrievelist),.progress = "text")
  names(results_as_lists) <- accessions

  if(!compact_results){
    results_as_lists
  }else{
    #Compact the result type dataframes into a single instance. Per accession counts in each column.
    all_results <- plyr::llply(retrievelist, function(y){
      tryCatch({
        r = lapply(results_as_lists, function(x){
          df <- as.data.frame(x[[y]])
          df
        })
        longform <- dplyr::bind_rows(r, .id = "analysis")
        cn <- colnames(longform)
        extras <- cn[!(cn %in% c("count","index_id", "analysis"))]
        final_df <- reshape2::dcast(longform, as.formula(paste(paste(extras,  collapse = " + "), " ~ analysis")), value.var = "count", fun.aggregate = sum)
        final_df}, error=function(x) NULL)
    })
  }
  names(all_results) <- retrievelist
  all_results
}




