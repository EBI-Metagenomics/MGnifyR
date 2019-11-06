# NOTES:

#Lack of consistency in the API - e.g. the only reason to search by study is to use the "centre_name" filter
# Only (afaics) "accession" accepts multiple csv values.




#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


library(httr)
library(phyloseq)
library(plyr)



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
analysis_filters = c('')

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
  new("mgnify_client", url=url, authtok=authtok, cachedir = cachepath)

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
#' @return \code{char} string for the API path of x/typeY
#' @examples
#' cl <- new("mgnify_client")
#' mgnify_get_x_for_y(cl, "MGYS00005126", "studies", "samples")
#'@export
mgnify_get_x_for_y <- function(client, x, typeX, typeY){
  #This one's easy - just rearrange the URLs
  if(typeX=="samples" & typeY %in% c("runs","studies")){
    paste( typeX,x,typeY, sep="/")
  }else if(typeX=="runs" & typeY == "analyses"){
    paste( typeX,x,typeY, sep="/")
  }
  else{
    #Do it the hard way with a callout
    json_dat = mgnify_query_json(client, paste(typeX, x, sep="/"))
    json_dat
    tgt_access = json_dat[[1]]$relationships[[typeY]]$data$id
    tgt_type = json_dat[[1]]$relationships[[typeY]]$data$type
    paste(tgt_type,tgt_access,sep="/")
    #substr(tgt_url, nchar(client@url) + 1, nchar(tgt_url))
  }
}


#Not exporting this - if people want to they can use the
# rjsonapi functionality
#'Coverting attribute lists to a single data.frame row
#'
#'\code{mgnify_attr_list_to_df} extracts the \code{attribute} entry in a JSONAPI result
#'@param json The \emph{raw} result list
#'@param metadata_key Optional extra key to parse subattributes from
#'@return data.frame containing a single row of metadata
#'
#'@export
mgnify_attr_list_to_df <- function (json, metadata_key=NULL ){

  attrlist=names(json$attributes)
  if (!is.null(metadata_key)){
    baseattrlist=attrlist[!(attrlist %in% c(metadata_key))]
    metaattrlist=json$attributes[[metadata_key]]
    metlist=sapply(metaattrlist, function(x) x$value)
    names(metlist)=sapply(metaattrlist, function(x) x$key)
    df = as.data.frame(t(unlist(c(json["attributes"][baseattrlist], metlist))), stringsAsFactors = F)
  }else{
    df = as.data.frame(t(unlist(json["attributes"])), stringsAsFactors = F)
  }
  df$accession <- json$id

  rownames(df)=df$accession
  df
}



#Internal function to actually perform the http request. Build up the URL then issues
#a GET, parsing the returned JSON into a nested list (uses \code{jsonlite} internally?)
#Previously cached results may be retrieved from disk without resorting to claling the MGnify server.

#'Low level MGnify API handler
#'
#'\code{mgnify_query_json} deals with handles the actual HTTP GET calls for the MGnifyR package, handling both pagination and local reuslt
#'caching. Although principally intended for internal MGnifyR use , it's exported for direct invocation.
#'
#'@param client MGnifyR client
#'@param path top level search point for the query. One of \code{biomes}, \code{samples}, \code{runs} etc.
#'@param qopts named list or vector containing options/filters to be URL encoded and appended to query as key/value pairs
#'@param maxhits Maxmium number of data entries to return. The actual number of hits returned may be higher than this value,
#'as this parameter only clamps after each full page is processed.
#'@param usecache Should successful queries be cached on disk locally? There are unresolved questions about whether this is
#'a sensible thing to do, but it remains as an option. It probably makes sense for single accession grabs, but not for
#'(filtered) queries - which are liable to change as new data is added to MGnify. Also caching only works for the first page.
#'@param Debug Should we print out lots of information while doing the grabbing?
#'@return \code{list} of results after pagination is dealt with.
#'@export
  mgnify_query_json <- function(client, path="biomes", qopts=NULL, maxhits=200, usecache = F, Debug=F){
  #Set up the base url
  fullurl = paste(client@url, path, sep="/")

  #convert to csv if filters are lists.
  #This doesn't check if they ~can~ be searched for in the API,
  #which is an issue since no error is returned by the JSON if the search
  #is invalid - we only get a
  tmpqopts = lapply(qopts,function(x) paste(x,collapse = ','))

  #Include the json and page position options
  full_qopts = as.list(c(format="json", tmpqopts, page=1))

  # Do we want to try and use a cache to speed things up?
  if(usecache){
    fname_list = c(path, names(unlist(full_qopts)), unlist(full_qopts))
    cache_fname = paste(fname_list,collapse = "_")
    cache_full_fname = paste(client@cache_dir, '/', cache_fname, '.RDS', sep="")
    if (file.exists(cache_full_fname)){
      data = readRDS(cache_full_fname)
    }else{
      res = httr::GET(url=fullurl, config(verbose=Debug), query=full_qopts )
      data <-httr::content(res)
      saveRDS(data, cache_full_fname)
    }
  }else{
    res = httr::GET(url=fullurl, config(verbose=Debug), query=full_qopts )
    data <-httr::content(res)
  }


  #Create something to store the returned data
  datlist=list()
  datlist[[1]] = data$data
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
      if(!is.null(maxhits)){
        curlen=sum(sapply(datlist, length))
        if (curlen > maxhits){
          break
        }
      }
    }
    unlist(datlist, recursive=F)
  }
  else{
    datlist
  }
}


#'#' Search MGnify database for studies, samples and runs
#'
#' \code{mgnify_query} is a flexible query function, harnessing the full power of the JSONAPI MGnify
#' search filters. Can be used for both metadata retrieval and
#' @param \code{mgnify_client} instance
#' @param \code{qtype} Type of objects to query. One of \code{studies},\code{samples},\code{runs} or
#' \code{analyses}
#' @param accession Either a single known MGnify accession identifier (of type \code{qtype}), or a list/vector
#' of accessions to query.
mgnify_query <- function(client, qtype="samples", accession=NULL, asDataFrame=F, ...){
  #Need to get around the lazy expansion in R in order to get a list
  a=accession
  arglist = as.list(match.call())[-1] # drop off the first entry, which is the name of the function
#  arglist

  arglist$accession=a

  #Filter the query options such that
  qopts = arglist[names(arglist) %in% query_filters[[qtype]]]
  non_qopts = arglist[!(names(arglist) %in% c(c("asDataFrame","qtype","client"),query_filters[[qtype]]))]

  cat(str(arglist))
  all_query_params = unlist(list(c(non_qopts,list(client=client, path=qtype, qopts=qopts))))

  cat(str(all_query_params))
  #Do the query
  #result = mgnify_query_json(client, path=qtype, qopts = qopts)
  result = do.call("mgnify_query_json", all_query_params)

  #Rename entries by accession
  id_list = lapply(result, function(x) x$id)
  names(result) = id_list

  if(asDataFrame){
    #Because metadata might not match across studies, the full dataframe is built by first building per-sample dataframes,
    # then using rbind.fill from plyr to combine. For ~most~ use cases the number of empty columns will hopefully
    # be minimal... because who's going to want cross study grabbing (?)
    dflist = lapply(result, function(r){
      df2 <- mgnify_attr_list_to_df(json = r, metadata_key = "sample-metadata")
      df2$biome = r$relationship$biome$data$id
      df2$study = r$relationship$studies$data$id
      df2$type = r$type
      rownames(df2)=df2$accession
      df2

    }
    )
    tryCatch(
      plyr::rbind.fill(dflist),
      error=function(e) dflist
    )
  }else{
    result
  }

}


#Using a previously retrieved (and possibly filtered) \code{mgnify_query} result(s), retrieve the corresponding
#analyses and
#associated dataset, for inclusion in a phyloseq object.


#mgnify_get_analyses(client, query_results)
#This does the heavy downloading of BIOM files for conversion into phyloseq.
#accessions can be either:
# - unnamed list or vector of run accessions:
# - named list returned from one of the mgnify_query_xxx functions
# - dataframe (rownames = accession IDs)
#
# pipeline_version is an extra filter to ensure only biomes of the same version get munged together
# otherwise by default the first pipeline version in the first sample/run will be used.
#


#' Retrieve \code{analyses} from MGnify and convert into \code{phyloseq} objects.
#'
#' \code{mgnify_get_runs_as_phyloseq} takes as input a previously determined data.frame or list from \code{mgnify_query}, and downloads all associated
#' analysis data (principally \code{.biom} OTU table output), before merging it with corresponding \code{sample}, \code{study} and \code{run} data into a
#' single \code{phyloseq} object. The \code{sample_data} in the resulting phyloseq object includes all \code{attributes} metadata from corresponding \code{
#' samples}, \code{studies} and \code{runs} entries
#'
#' @param client MGnifyR client
#' @param accessions data.frame containing (at least) two columns: \code{accession} - MGnify accession id (sample, run, study etc), and
#' \code{type} - type of accession in \code{accession} column. Lists of accessions will be acceptable in the future. \code{type} defaults to \code{samples}
#' if type column is absent. In most cases, the \code{accession} input will be the output of a previous \code{mgnify_query} call.
#' @param downloadDIR Location to store retrieved \code{.biom} files. Can be resued between session to reduce load and speed up analysis.
#' @param use_downloads Try to use previously retrieved \code{.biom} files, or start from scratch.
#' @param pipeline_version Unimplemented - filtering ~should~ be done at the data.frame level, before we get to this point?
#' @param usecache Should JSON API queries be cached locally for performance?
#' @return \code{phyloseq} object with filled \code{sample_data} and \code{otu_table} slots.
#' @example
#'
#'@export
mgnify_get_runs_as_phyloseq <- function(client=NULL, accessions=NULL, downloadDIR='./tmpdownloads', use_downloads=T, pipeline_version, usecache=T ){
  #These must be RUN accessions
  dir.create(downloadDIR, showWarnings = F)
  grabtype="unk"
  #At the end of this, we want a list of run accessions to grab # unimplemented for now... just data.frame works atm
  if (class(accessions) == "list"){
    #Unnamed list
    if(is.null(names(accessions))){
      grabtype="raw_run"
    }
    else{
      grabtype="namedlist"
    }

  }else if(class(accessions) =="data.frame"){
    #Same as above, this time using the rownames as accession numbers, and the "type" column to figure out where to go
    grabtype="data.frame"
    full_ps_list <- list()
    for (cur_line in seq(nrow(accessions))){
      if ("type" %in% colnames(accessions)){
        curtype=accessions[cur_line, "type"]
        curaccess=accessions$accession[cur_line]
      }else{
        #Assume that they're samples
        curtype = "samples"
        curaccess = accessions$accession[cur_line]
      }
      #We need to get to the "runs" section, and then to the analysis:
      #Each "run" has only one "analyses", so makes sense to go through "runs"
      #We might already be there of course...
      if (curtype != "run"){
        runpath=mgnify_get_x_for_y(client, curaccess, curtype, "runs" )
      } else{
        runpath=paste("runs",curaccess)
      }


      #Retrieve the run data
      run_data <- mgnify_query_json(client, runpath, usecache = usecache)
      #get the sample and study data as well. This repeats some earlier calls, but makes it easier
      #to code up:
      sample_data <- mgnify_query_json(client, paste("samples",run_data[[1]]$relationships$sample$data$id, sep="/"),usecache = usecache)
      study_data <- mgnify_query_json(client, paste("studies",run_data[[1]]$relationships$study$data$id, sep="/"),usecache = usecache)

      samp_attr_df <- mgnify_attr_list_to_df(sample_data[[1]], "sample-metadata")
      study_attr_df <- mgnify_attr_list_to_df(study_data[[1]])

      #Depending how we got there, run_data might be length > 1, so a quick rename by accession should make
      #it easier later
      names(run_data) <- sapply(run_data, `[`, "id")

      #Each run will have an analysis entry:
      analyses_phyloseqs <- sapply(names(run_data), function(r){
        analyse_path = paste("runs",run_data[[r]]$id,"analyses", sep="/")
        cat(str(analyse_path))
        analysis_data <- mgnify_query_json(client, analyse_path,usecache = usecache)
        #Build up a dataframe of attributes
        analysis_attr_df <- mgnify_attr_list_to_df(analysis_data[[1]], metadata_key = "analysis-summary")
        analysis_attr_df

        #Get the download page json
        analysis_downloads <- mgnify_query_json(cl, paste("analyses", analysis_attr_df$accession, "downloads", sep="/"),usecache = usecache)
        #find out where our biom file is:
        biom_url <- analysis_downloads[grepl('JSON Biom', sapply(analysis_downloads, function(x) x$attributes$`file-format`$name))][[1]]$links$self
        fname=tail(strsplit(biom_url, '/')[[1]], n=1)
        biom_path = paste(downloadDIR, fname, sep="/")
        if (! file.exists(biom_path) | !use_downloads ){
          httr::GET(biom_url, write_disk(biom_path, overwrite = T))
        }
        #Load in the phlyloseq object
        psobj <- phyloseq::import_biom(biom_path)
        #The biom files have a single column of "sa1". It's rewritten as sample_run_analysis accession, with
        # the original value stored in the sample_data (just in case it changes between pipelines)
        orig_samp_name <- sample_names(psobj)[[1]]
        newsampname <- paste(analysis_attr_df$accession)
        sample_names(psobj) <- newsampname

        #saveRDS(psobj, paste(biom_path,".RDS",sep=""))
        colnames(samp_attr_df) <- paste("SAMPLE_",colnames(samp_attr_df),sep='')
        colnames(study_attr_df) <- paste("STUDY_",colnames(study_attr_df),sep='')
        colnames(analysis_attr_df) <- paste("ANALYSIS_",colnames(analysis_attr_df),sep='')
        full_samp_data <- cbind(samp_attr_df, study_attr_df, analysis_attr_df)
        rownames(full_samp_data) <- newsampname
        sample_data(psobj) <- full_samp_data
        #saveRDS(psobj, paste(biom_path,".RDS2",sep=""))
        psobj
      })

      full_ps_list[[cur_line]] = analyses_phyloseqs
    }
  }
  names(full_ps_list) <- NULL
  full_ps <- do.call(phyloseq::merge_phyloseq, unlist(full_ps_list))
  #For some reason merge_phyloseq messes up sample data, so we need to rebuild the data.frame and
  #add it to the phyloseq object:
  sampdatlist <- lapply(unlist(full_ps_list), function(x) as.data.frame(phyloseq::sample_data(x)))
  full_sample_df <-  do.call(plyr::rbind, sampdatlist)
  rownames(full_sample_df) <- full_sample_df$ANALYSIS_accession
  phyloseq::sample_data(full_ps) <- full_sample_df
  full_ps
}



