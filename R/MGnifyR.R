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
    samplist = sapply(result, function(r){
      attrlist=names(r$attributes)
      baseattrlist=attrlist[!(attrlist %in% "sample-metadata")]
      metaattrlist=r$attributes$`sample-metadata`
      metlist=sapply(metaattrlist, function(x) x$value)
      names(metlist)=sapply(metaattrlist, function(x) x$key)
      df = as.data.frame(t(unlist(c(r$attributes[baseattrlist], metlist))), stringsAsFactors = F)
      df$biome = r$relationship$biome$data$id
      df$study = r$relationship$studies$data[[1]]$id
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


#Retrieves runs/analyses
mgnify_query_studies<- function(client, accession=NULL, asDataFrame=F, ...){
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
    samplist = sapply(result, function(r){
      attrlist=names(r$attributes)
      baseattrlist=attrlist[!(attrlist %in% "sample-metadata")]
      metaattrlist=r$attributes$`sample-metadata`
      metlist=sapply(metaattrlist, function(x) x$value)
      names(metlist)=sapply(metaattrlist, function(x) x$key)
      df = as.data.frame(t(unlist(c(r$attributes[baseattrlist], metlist))), stringsAsFactors = F)
      df$biome = r$relationship$biome$data$id
      df$study = r$relationship$studies$data[[1]]$id
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





#res=mgnify_query_samples(cl, biome_name="root:Engineered:Wastewater", asDataFrame = T)#

#res2=mgnify_query_samples(cl, accession = res$accession, asDataFrame = T)

#res$accession




