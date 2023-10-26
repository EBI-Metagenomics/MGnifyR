# All generic methods are listed here

#' @export
setGeneric("url", signature = c("x"), function(x) standardGeneric("url"))

#' @export
setGeneric("authTok", signature = c("x"), function(x)
    standardGeneric("authTok"))

#' @export
setGeneric("cacheDir", signature = c("x"), function(x)
    standardGeneric("cacheDir"))

#' @export
setGeneric("warnings", signature = c("x"), function(x)
    standardGeneric("warnings"))

#' @export
setGeneric("useMemCache", signature = c("x"), function(x)
    standardGeneric("useMemCache"))

#' @export
setGeneric("memCache", signature = c("x"), function(x)
    standardGeneric("memCache"))

#' @export
setGeneric("clearCache", signature = c("x"), function(x)
    standardGeneric("clearCache"))

#' @export
setGeneric(
    "url<-", signature = c("x"), function(x, value) standardGeneric("url<-"))

#' @export
setGeneric("authTok<-", signature = c("x"), function(x, value)
    standardGeneric("authTok<-"))

#' @export
setGeneric("cacheDir<-", signature = c("x"), function(x, value)
    standardGeneric("cacheDir<-"))

#' @export
setGeneric("warnings<-", signature = c("x"), function(x, value)
    standardGeneric("warnings<-"))

#' @export
setGeneric("useMemCache<-", signature = c("x"), function(x, value)
    standardGeneric("useMemCache<-"))

#' @export
setGeneric("memCache<-", signature = c("x"), function(x, value)
    standardGeneric("memCache<-"))

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
