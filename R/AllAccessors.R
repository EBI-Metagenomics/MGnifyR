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
#' databaseUrl(mg)
#' showWarnings(mg) <- FALSE
#'
#' @name MgnifyClient-accessors
NULL

#' @rdname MgnifyClient-accessors
#' @include AllClasses.R AllGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "databaseUrl", signature = c(x = "MgnifyClient"),
    function(x){ x@databaseUrl })

#' @rdname MgnifyClient-accessors
#' @include AllClasses.R AllGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "authTok", signature = c(x = "MgnifyClient"), function(x){ x@authTok })

#' @rdname MgnifyClient-accessors
#' @include AllClasses.R AllGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "useCache", signature = c(x = "MgnifyClient"),
    function(x){ x@useCache })

#' @rdname MgnifyClient-accessors
#' @include AllClasses.R AllGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "cacheDir", signature = c(x = "MgnifyClient"), function(x){ x@cacheDir })

#' @rdname MgnifyClient-accessors
#' @include AllClasses.R AllGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "showWarnings", signature = c(x = "MgnifyClient"),
    function(x){ x@showWarnings })

#' @rdname MgnifyClient-accessors
#' @include AllClasses.R AllGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "useMemCache", signature = c(x = "MgnifyClient"),
    function(x){ x@useMemCache })

#' @rdname MgnifyClient-accessors
#' @include AllClasses.R AllGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "memCache", signature = c(x = "MgnifyClient"), function(x){ x@memCache })

#' @rdname MgnifyClient-accessors
#' @include AllClasses.R AllGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "clearCache", signature = c(x = "MgnifyClient"),
    function(x){ x@clearCache })

#' @rdname MgnifyClient-accessors
#' @include AllClasses.R AllGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "verbose", signature = c(x = "MgnifyClient"),
    function(x){ x@verbose })

#' @rdname MgnifyClient-accessors
#' @include AllClasses.R AllGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "databaseUrl<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, databaseUrl = value) })

#' @rdname MgnifyClient-accessors
#' @include AllClasses.R AllGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "authTok<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, authTok = value) })

#' @rdname MgnifyClient-accessors
#' @include AllClasses.R AllGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "useCache<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, useCache = value) })

#' @rdname MgnifyClient-accessors
#' @include AllClasses.R AllGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "cacheDir<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, cacheDir = value) })

#' @rdname MgnifyClient-accessors
#' @include AllClasses.R AllGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "showWarnings<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, showWarnings = value) })

#' @rdname MgnifyClient-accessors
#' @include AllClasses.R AllGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "useMemCache<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, useMemCache = value) })

#' @rdname MgnifyClient-accessors
#' @include AllClasses.R AllGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "memCache<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, memCache = value) })

#' @rdname MgnifyClient-accessors
#' @include AllClasses.R AllGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "clearCache<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, clearCache = value) })

#' @rdname MgnifyClient-accessors
#' @include AllClasses.R AllGenerics.R MgnifyClient.R utils.R
#' @export
setMethod(
    "verbose<-", signature = c(x = "MgnifyClient"),
    function(x, value){ BiocGenerics:::replaceSlots(x, verbose = value) })
