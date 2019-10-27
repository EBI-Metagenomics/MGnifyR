# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


require(httr)
require(phyloseq)
require(ape)



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
              'metadata_key','metadata_value_gte','metadata_value_lte','metadata_value','sample_accession','study_accession',
              'include')


#Base class for retrieving URL queries from MGnify
#   support for paging

#Functions:
#Generic function to retrieve (paginated) data from mgnify API
#Classes

setClass("mgnify_query_getter", slots=list(url = "character", qopts = "list", maxreturn="numeric" ))
setClass("mgnify_search_class", slots=list(baseurl = "character", qfilters = "list", qfiltersvalid="list", qtype="character" ))



mgnify_get_json_data <- function(object,...){}

#generic function to build up the query array, parse the output ad return a dataframe
mgnify_search <- function(object){
  #Set up the base url
  url = paste(object@baseurl, object@qtype, sep="/")
  cat(url)
  qopts = qopts=c(format="json", object@qfilters)


  m <- new("mgnify_query_getter", url=url, qopts=qopts)

  dat <- mgnify_get_json_data(m, maxreturn=200)
  dat
  # build up the query parameter string
  # qstr_list = lapply(names(object@qfilters), function(n){paste(n,object@qfilters[[n]], sep="=")})
  #result = as.data.frame(t(as.data.frame(sapply(dat, function(x) {unlist(c(id=x[["id"]], x[["attributes"]]))}))))
  #rownames(result) <- dat$id
}


#Query a mgnify endpoint and return a list of $data entries - handles pagination
setMethod("mgnify_get_json_data","mgnify_query_getter", function(object, maxreturn=NULL){
  qopts = c(object@qopts, page=1)
  #cat(str(qopts))
  # Do a first grab of the data
  data <- content(GET(object@url, config(verbose=T), query=qopts ))

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
      qopts$page=p
      curd = content(GET(object@url, config(verbose=T), query=qopts ))
      datlist[[p]] = curd$data
      #Check to see if we've pulled enough entries
      if(!is.null(maxreturn)){
        curlen=sum(sapply(datlist, length))
        if (curlen > maxreturn){
          break
        }
      }
    }
    unlist(datlist, recursive=F)
  }
  else{
    datlist
  }

})


searchSamples


newsearch<-new("mgnify_search_class", baseurl=baseurl, qtype="samples", qfilters=list(experiment_type="amplicon"))
res = mgnify_search(newsearch)





res

