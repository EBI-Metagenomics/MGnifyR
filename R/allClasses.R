# MgnifyClient class and its accessors

#' A MgnifyClient object
#'
#' @details An object that are required by functions of MGnifyR package.
#'
#' @slot databaseUrl A single character value specifying an URL address of
#' database.
#'
#' @slot authTok A single character value specifying authentication token.
#' 
#' @slot useCache A single boolean value specifying whether to use cache.
#'
#' @slot cacheDir A single character value specifying cache directory.
#'
#' @slot showWarnings A single boolean value specifying whether to show
#' warnings.
#'
#' @slot useMemCache A single boolean value specifying whether to use on-disk
#' memory.
#'
#' @slot memCache A single character value specifying on-disk memory directory.
#'
#' @slot clearCache A single boolean value specifying whether to clear cache.
#' 
#' @slot verbose A single boolean value specifying whether to show messages.
#'
#' @section Constructor:
#' See  \code{\link{MgnifyClient}} for constructor.
#' 
#' @section Accessor:
#' See \code{\link{MgnifyClient-accessors}} for accessor functions.
#'
#' @name MgnifyClient
NULL

#' @rdname MgnifyClient
#' @importFrom httr POST
#' @importFrom httr content
#' @exportClass MgnifyClient
setClass(
    "MgnifyClient", representation(
        databaseUrl = "character",
        authTok = "character",
        useCache = "logical",
        cacheDir = "character",
        showWarnings = "logical",
        useMemCache = "logical",
        memCache = "list",
        clearCache = "logical",
        verbose = "logical"),
    prototype = list(
        databaseUrl = "https://www.ebi.ac.uk/metagenomics/api/v1",
        authTok = NULL,
        useCache = FALSE,
        cacheDir = NULL,
        useMemCache = FALSE,
        memCache = list(),
        clearCache = FALSE,
        verbose = TRUE))
