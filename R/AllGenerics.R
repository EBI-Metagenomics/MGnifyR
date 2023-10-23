# All generic methods are listed here

#' @rdname doQuery
#' @include MgnifyClient.R
#' @export
setGeneric(
    "doQuery", signature = c("x"), function(x, ...) standardGeneric("doQuery"))

#' @rdname getFile
#' 
#' @include MgnifyClient.R
#' @export
setGeneric(
    "getFile", signature = c("x"), function(x, ...) standardGeneric("getFile"))

#' @rdname getFile
#' @include MgnifyClient.R
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
