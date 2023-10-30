#' MgnifyClient accessors and mutators
#'
#' @details
#' These functions are for fetching and mutating slots of
#' \code{MgnifyClient} object.
#'
#' @param x A \code{MgnifyClient} object.
#'
#' @param value A value to be added to a certain slot.
#'
#' @return A value of MgnifyClient object or nothing.
#'
#' @examples
#' mg <- MgnifyClient()
#' 
#' url(mg)
#' warnings(mg) <- FALSE
#'
#' @name MgnifyClient-accessors
NULL

#' @rdname MgnifyClient-accessors
#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod("url", signature = c(x = "MgnifyClient"), function(x){ x@url })

#' @rdname MgnifyClient-accessors
#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "authTok", signature = c(x = "MgnifyClient"), function(x){ x@authTok })

#' @rdname MgnifyClient-accessors
#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "cacheDir", signature = c(x = "MgnifyClient"), function(x){ x@cacheDir })

#' @rdname MgnifyClient-accessors
#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "warnings", signature = c(x = "MgnifyClient"), function(x){ x@warnings })

#' @rdname MgnifyClient-accessors
#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "useMemCache", signature = c(x = "MgnifyClient"),
    function(x){ x@useMemCache })

#' @rdname MgnifyClient-accessors
#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "memCache", signature = c(x = "MgnifyClient"), function(x){ x@memCache })

#' @rdname MgnifyClient-accessors
#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "clearCache", signature = c(x = "MgnifyClient"),
    function(x){ x@clearCache })

#' @rdname MgnifyClient-accessors
#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "url<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, url = value) })

#' @rdname MgnifyClient-accessors
#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "authTok<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, authTok = value) })

#' @rdname MgnifyClient-accessors
#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "cacheDir<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, cacheDir = value) })

#' @rdname MgnifyClient-accessors
#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "warnings<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, warnings = value) })

#' @rdname MgnifyClient-accessors
#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "useMemCache<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, useMemCache = value) })

#' @rdname MgnifyClient-accessors
#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "memCache<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, memCache = value) })

#' @rdname MgnifyClient-accessors
#' @include allClass.R allGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "clearCache<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, clearCache = value) })
