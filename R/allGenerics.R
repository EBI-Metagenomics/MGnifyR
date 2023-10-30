# All generic methods are listed here

#' @rdname MgnifyClient-accessors
#' @export
setGeneric("url", signature = c("x"), function(x) standardGeneric("url"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric("authTok", signature = c("x"), function(x)
    standardGeneric("authTok"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric("cacheDir", signature = c("x"), function(x)
    standardGeneric("cacheDir"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric("warnings", signature = c("x"), function(x)
    standardGeneric("warnings"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric("useMemCache", signature = c("x"), function(x)
    standardGeneric("useMemCache"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric("memCache", signature = c("x"), function(x)
    standardGeneric("memCache"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric("clearCache", signature = c("x"), function(x)
    standardGeneric("clearCache"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric(
    "url<-", signature = c("x"), function(x, value) standardGeneric("url<-"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric("authTok<-", signature = c("x"), function(x, value)
    standardGeneric("authTok<-"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric("cacheDir<-", signature = c("x"), function(x, value)
    standardGeneric("cacheDir<-"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric("warnings<-", signature = c("x"), function(x, value)
    standardGeneric("warnings<-"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric("useMemCache<-", signature = c("x"), function(x, value)
    standardGeneric("useMemCache<-"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric("memCache<-", signature = c("x"), function(x, value)
    standardGeneric("memCache<-"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric("clearCache<-", signature = c("x"), function(x, value)
    standardGeneric("clearCache<-"))

#' @rdname doQuery
#' @export
setGeneric(
    "doQuery", signature = c("x"), function(x, ...) standardGeneric("doQuery"))

#' @rdname getFile
#' @export
setGeneric(
    "getFile", signature = c("x"), function(x, ...) standardGeneric("getFile"))

#' @rdname getFile
#' @export
setGeneric(
    "searchFile", signature = c("x"), function(x, ...)
        standardGeneric("searchFile"))

#' @rdname getMetadata
#' @export
setGeneric(
    "getMetadata", signature = c("x"), function(x, ...)
        standardGeneric("getMetadata"))

#' @rdname getResult
#' @export
setGeneric(
    "getResult", signature = c("x"), function(x, ...)
        standardGeneric("getResult"))

#' @rdname searchAnalysis
#' @export
setGeneric(
    "searchAnalysis", signature = c("x"), function(x, ...)
        standardGeneric("searchAnalysis"))
