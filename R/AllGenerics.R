# All generic methods are listed here

#' @rdname MgnifyClient-accessors
#' @export
setGeneric(
    "databaseUrl", signature = c("x"), function(x)
        standardGeneric("databaseUrl"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric(
    "authTok", signature = c("x"), function(x) standardGeneric("authTok"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric(
    "useCache", signature = c("x"), function(x) standardGeneric("useCache"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric(
    "cacheDir", signature = c("x"), function(x) standardGeneric("cacheDir"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric(
    "showWarnings", signature = c("x"), function(x)
        standardGeneric("showWarnings"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric(
    "clearCache", signature = c("x"), function(x) standardGeneric("clearCache"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric(
    "verbose", signature = c("x"), function(x) standardGeneric("verbose"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric(
    "databaseUrl<-", signature = c("x"), function(x, value)
        standardGeneric("databaseUrl<-"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric(
    "authTok<-", signature = c("x"), function(x, value)
        standardGeneric("authTok<-"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric(
    "useCache<-", signature = c("x"), function(x, value)
        standardGeneric("useCache<-"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric(
    "cacheDir<-", signature = c("x"), function(x, value)
        standardGeneric("cacheDir<-"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric(
    "showWarnings<-", signature = c("x"), function(x, value)
        standardGeneric("showWarnings<-"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric(
    "clearCache<-", signature = c("x"), function(x, value)
        standardGeneric("clearCache<-"))

#' @rdname MgnifyClient-accessors
#' @export
setGeneric(
    "verbose<-", signature = c("x"), function(x, value)
        standardGeneric("verbose<-"))

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
