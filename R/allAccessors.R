# Here are listed all MgnifyClient class accessors and mutators.

################################## ACCESSORS ###################################

#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod("url", signature = c(x = "MgnifyClient"), function(x){ x@url })

#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "authTok", signature = c(x = "MgnifyClient"), function(x){ x@authTok })

#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "cacheDir", signature = c(x = "MgnifyClient"), function(x){ x@cacheDir })

#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "warnings", signature = c(x = "MgnifyClient"), function(x){ x@warnings })

#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "useMemCache", signature = c(x = "MgnifyClient"),
    function(x){ x@useMemCache })

#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "memCache", signature = c(x = "MgnifyClient"), function(x){ x@memCache })

#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "clearCache", signature = c(x = "MgnifyClient"),
    function(x){ x@clearCache })

################################### MUTATORS ###################################

#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "url<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, url = value) })

#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "authTok<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, authTok = value) })

#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "cacheDir<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, cacheDir = value) })

#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "warnings<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, warnings = value) })

#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "useMemCache<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, useMemCache = value) })

#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "memCache<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, memCache = value) })

#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "clearCache<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, clearCache = value) })
