# NOTES:

#Lack of consistency in the API - e.g. the only reason to search by study is to use the "centre_name" filter
# Only (afaics) "accession" accepts multiple csv values.




#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#library(httr)
require(httr)
require(phyloseq)
require(ape)
require(plyr)



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
##'Coverting attribute lists to a single data.frame row
##'
##'\code{mgnify_attr_list_to_df} extracts the \code{attribute} entry in a JSONAPI result
##'@param json The \emph{raw} result list
##'@param metadata_key Optional extra key to parse child entries and include in the output
##'@return 1xn data.frame
##'@export
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



mgnify_query_json <- function(client, path="biomes", qopts=NULL, maxhits=200, ...){
  #Set up the base url
  fullurl = paste(client@url, path, sep="/")
  #cat(fullurl)

  #convert to csv if filters are lists.
  #This doesn't check if they ~can~ be
  tmpqopts = lapply(qopts,function(x) paste(x,collapse = ','))

  #Include the json and page position options
  full_qopts = as.list(c(format="json", tmpqopts, page=1))
  cat(str(full_qopts))

  # Do a first grab of the data
  res = GET(url=fullurl, config(verbose=T), query=full_qopts )
  cat(str(res))
  data <-content(res)

  #Create something to store the returned data
  datlist=list()
  datlist[[1]] = data$data
  # Check to see if there's pagination required
  if ("meta" %in% names(data)){
    #Yes, paginate
    pstart = as.numeric(data$meta$pagination$page)
    pend   = as.numeric(data$meta$pagination$pages)

    for (p in seq(pstart+1,pend)){  # We've already got the first one
      cat(p)
      #curd <- rjson::fromJSON(file=paste(object@url, "&page=", p, sep=""))
      full_qopts$page=p
      curd = content(GET(fullurl, config(verbose=T), query=full_qopts ))
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



#'@export
mgnify_query <- function(client, qtype="samples", accession=NULL, asDataFrame=F, ...){
  #Need to get around the lazy expansion in R in order to get a list
  a=accession
  arglist =as.list(match.call())
  arglist$accession=a

  #Filter the query options such that
  qopts = arglist[names(arglist) %in% query_filters[[qtype]]]

  #Do the query
  result = mgnify_query_json(client, path=qtype, qopts = qopts)

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



#Search by "sample":
#various parameters to search by


#Function to retrieve sample information as a data.frame.
#Two forms:
# - provide a (list) of sample accession numbers, possible from a previous search (e.g. by project)
# - populate a named list of filters (see sample_filters above
# R could really do with dictionaries because this'd be a lot easier then.
# By default, because samples may be associated with multiple studies, and may
# be analysed during multiple "runs", this'll return a named list of lists. Top-level
# names are the sample accession, with level two entries corresponding to individual attributes.
# Because most samples are only belonging to one project, and only get analysed once, the "asDataFrame"
# option allows forced coercion into a data.frame by selecting the first entry in any lists. This should work
# for most use cases.
# Attributes and metadata_attributes are expanded as columns in the data.frames.

#Annoyingly this won't work with a list of project accessions (which IMHO it should). thus you need to iterate over a
#project list externally to get it to work. Which then makes joining the results together a bit awkward. Hey ho.
#' @export
mgnify_query_samples <-
  function(client, accession=NULL, asDataFrame=F, ...){
  #Is there a proper way to do this? F'in lazy evaluation:
  a=accession
  arglist =as.list(match.call())
  arglist$accession=a
  #cat(str(names(arglist)))
  #COnvert the extra arguments into valid MGnify query filter arguments
  qopts = arglist[names(arglist) %in% sample_filters]


  result = mgnify_query_json(client, path="samples", qopts = qopts)

  samp_id_list = lapply(result, function(x) x$id)
  names(result) = samp_id_list

  if(asDataFrame){
    #Because metadata might not match across studies, the full dataframe is built by first building per-sample dataframes,
    # then using rbind.fill from plyr to combine. For ~most~ use cases the number of empty columns will hopefully
    # be minimal... because who's going to want cross study grabbing (?)
    samplist = lapply(result, function(r){
      df2 <- mgnify_attr_list_to_df(json = r, metadata_key = "sample-metadata")
      df2$biome = r$relationship$biome$data$id
      df2$study = r$relationship$studies$data$id
      df2$type = r$type
      rownames(df2)=df2$accession
      df2

    }
    )
    tryCatch(
      plyr::rbind.fill(samplist),
      error=function(e) samplist
    )
  }else{
    result
  }

}


#Retrieves studies - again, only study accessions can be lists.

mgnify_query_studies<- function(client, accession=NULL, asDataFrame=F, ...){
  #Is there a proper way to do this? F'in lazy evaluation:
  a=accession
  arglist =as.list(match.call())
  arglist$accession=a
  #cat(str(names(arglist)))
  #Convert the extra arguments into valid MGnify query filter arguments
  qopts = arglist[names(arglist) %in% study_filters]


  result = mgnify_query_json(client, path="studies", qopts = qopts)

  study_id_list = lapply(result, function(x) x$id)
  names(result) = study_id_list

  if(asDataFrame){
    #Studies don't have metadata, so we're just returning the list of attributes.
    studylist = lapply(result, function(r){
      attrlist=names(r$attributes)
      baseattrlist=attrlist

      df = as.data.frame(t(unlist(r$attributes[baseattrlist])), stringsAsFactors = F)
      df$biome = r$relationship$biome$data$id
      df$study = r$relationship$studies$data[[1]]$id
      rownames(df)=df$accession
      df$type = r$type
      df

    }
    )
    tryCatch(
      plyr::rbind.fill(studylist),
      error=function(e) studylist
    )
  }else{
    result
  }

}

#Does essentially the same thing as the two other functions above, but for runs instead of
#Or at least it will do once it's implemented...

mgnify_query_runs<- function(client, accession=NULL, asDataFrame=F, ...){}




#This does the heavy downloading of BIOM files for conversion into phyloseq.
#accessions can be either:
# - unnamed list or vector of run accessions:
# - named list returned from one of the mgnify_query_xxx functions
# - dataframe (rownames = accession IDs)
#
# pipeline_version is an extra filter to ensure only biomes of the same version get munged together
# otherwise by default the first pipeline version in the first sample/run will be used.
#
mgnify_get_runs_as_phyloseq <- function(client=NULL, accessions=NULL, downloadDIR='./tmpdownloads', use_downloads=T, pipeline_version){
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

      cat(runpath)
      #Retrieve the run data
      run_data <- mgnify_query_json(client, runpath)
      #get the sample and study data as well. This repeats some earlier calls, but makes it easier
      #to code up:
      sample_data <- mgnify_query_json(client, paste("samples",run_data[[1]]$relationships$sample$data$id, sep="/"))
      study_data <- mgnify_query_json(client, paste("studies",run_data[[1]]$relationships$study$data$id, sep="/"))

      samp_attr_df <- mgnify_attr_list_to_df(sample_data[[1]], "sample-metadata")
      study_attr_df <- mgnify_attr_list_to_df(study_data[[1]])

      #Depending how we got there, run_data might be length > 1, so a quick rename by accession should make
      #it easier later
      names(run_data) <- sapply(run_data, `[`, "id")
      cat(str(run_data))
      #Each run will have an analysis entry:
      analyses_phyloseqs <- sapply(names(run_data), function(r){
        analyse_path = paste("runs",run_data[[r]]$id,"analyses", sep="/")
        cat(str(analyse_path))
        analysis_data <- mgnify_query_json(client, analyse_path)
        #Build up a dataframe of attributes
        analysis_attr_df <- mgnify_attr_list_to_df(analysis_data[[1]], metadata_key = "analysis-summary")
        analysis_attr_df

        #Get the download page json
        analysis_downloads <- mgnify_query_json(cl,
                                                paste("analyses", analysis_attr_df$accession, "downloads", sep="/"))
        #find out where our biom file is:
        biom_url <- analysis_downloads[grepl('JSON Biom', sapply(analysis_downloads, function(x) x$attributes$`file-format`$name))][[1]]$links$self
        fname=tail(strsplit(biom_url, '/')[[1]], n=1)
        biom_path = paste(downloadDIR, fname, sep="/")
        if (! file.exists(biom_path) | !use_downloads ){
          GET(biom_url, write_disk(biom_path, overwrite = T))
        }
        #Load in the phlyloseq object
        psobj <- import_biom(biom_path)
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
  full_ps <- do.call(merge_phyloseq, unlist(full_ps_list))
  #For some reason merge_phyloseq messes up sample data, so we need to rebuild the data.frame and
  #add it to the phyloseq object:
  sampdatlist <- lapply(unlist(full_ps_list), function(x) as.data.frame(sample_data(x)))
  full_sample_df <-  do.call(rbind, sampdatlist)
  rownames(full_sample_df) <- full_sample_df$ANALYSIS_accession
  sample_data(full_ps) <- full_sample_df
  full_ps
}


#sample_df <- mgnify_query_samples(cl, study_accession = "MGYS00005120", asDataFrame = T)
#mgnify_get_runs_as_phyloseq(cl, accessions = sample_df)


#}

