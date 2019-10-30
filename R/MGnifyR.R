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


#Base class for retrieving URL queries from MGnify
#   support for paging

#Functions:
#Generic function to retrieve (paginated) data from mgnify API
#Classes

setClass("mgnify_client",
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

#Generic function to run a MGNify query and return JSON parsed as a list
#Requires mgnify_client instance as an argument, with optional parameters for
#target path, as well as query options (as a named list)
#Mostly used internally by more specific search code
mgnify_query_json <- function(client, path="biomes", qopts=NULL, maxhits=200, ...){
  #Set up the base url
  fullurl = paste(client@url, path, sep="/")
  #cat(fullurl)

  #convert to csv if filters are lists.
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
mgnify_query_samples<- function(client, accession=NULL, asDataFrame=F, ...){
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
      attrlist=names(r$attributes)
      baseattrlist=attrlist[!(attrlist %in% "sample-metadata")]
      metaattrlist=r$attributes$`sample-metadata`
      metlist=sapply(metaattrlist, function(x) x$value)
      names(metlist)=sapply(metaattrlist, function(x) x$key)
      df = as.data.frame(t(unlist(c(r$attributes[baseattrlist], metlist))), stringsAsFactors = F)
      df$biome = r$relationship$biome$data$id
      df$study = r$relationship$studies$data[[1]]$id
      df$type = r$type
      rownames(df)=df$accession
      df

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
mgnify_query_runs<- function(client, accession=NULL, asDataFrame=F, ...){}


#helper function for getting relative paths in the API
#Not everything is implemented here - just what we
#need to get to the download or run areas
#Given an accession x, we want to get the link to get the url for the
#coresponding typeY

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

#This does the heavy downloading of BIOM files for conversion into phyloseq.
#accessions can be either:
# - unnamed list or vector of run accessions:
# - named list returned from one of the mgnify_query_xxx functions
# - dataframe (rownames = accession IDs)
#
# pipeline_version is an extra filter to ensure only biomes of the same version get munged together
# otherwise by default the first pipeline version in the first sample/run will be used.
#
mgnify_get_runs_as_phyloseq <- function(client=NULL, accessions=NULL, downloadDIR, pipeline_version){
  #These must be RUN accessions
  grabtype="unk"
  #At the end of this, we want a list of run accessions to grab
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
    for (cur_line in nrow(accessions)){
      if ("type" %in% colnames(accessions)){
        curtype=accessions[cur_line, "type"]
        curaccess=accessions$accession[cur_line]
      }else{
        #Assume that they're samples
        curtype = "sample"
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
      run_data <- mgnify_query_json(client, runpath)
      #Depending how we got there, run_data might be length > 1, so a quick rename by accession should make
      #it easier later
      names(run_data) <- sapply(run_data, `[`, "id")
    cat(str(run_data))
      #Each run will have an analysis entry:
      analyses_accessions <- sapply(names(run_data), function(r){
        analyse_path = paste("runs",run_data[[r]]$id,"analyses", sep="/")
        cat(str(analyse_path))
        analysis_data <- mgnify_get_json_data(client, analyse_path)
        analysis_data$id
      })
    }

  }
  analyses_accessions

}


#sample_df <- mgnify_query_samples(cl, study_accession = "MGYS00005120", asDataFrame = T)
#mgnify_get_runs_as_phyloseq(cl, accessions = sample_df)


#}

