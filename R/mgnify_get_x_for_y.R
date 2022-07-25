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
## @export
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
