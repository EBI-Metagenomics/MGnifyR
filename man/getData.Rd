% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/getData.R
\name{getData}
\alias{getData}
\alias{getData,MgnifyClient-method}
\title{Versatile function to retrieve raw results}
\usage{
getData(x, ...)

\S4method{getData}{MgnifyClient}(x, type, accession.type = NULL, accession = NULL, as.df = TRUE, ...)
}
\arguments{
\item{x}{A \code{MgnifyClient} object.}

\item{...}{optional arguments fed to internal functions.}

\item{type}{A single character value specifying the type of data retrieve.
Must be one of the following options: \code{studies}, \code{samples},
\code{runs}, \code{analyses}, \code{biomes}, \code{assemblies},
\code{super-studies}, \code{experiment-types}, \code{pipelines},
\code{pipeline-tools}, \code{publications}, \code{genomes},
\code{genome-search}, \code{genome-search/gather}, \code{genome-catalogues},
\code{genomeset}, \code{cogs}, \code{kegg-modules}, \code{kegg-classes},
\code{antismash-geneclusters}, \code{annotations/go-terms},
\code{annotations/interpro-identifiers}, \code{annotations/kegg-modules},
\code{annotations/pfam-entries}, \code{annotations/kegg-orthologs},
\code{annotations/genome-properties},
\code{annotations/antismash-gene-clusters}, \code{annotations/organisms}, or
\code{mydata}.}

\item{accession.type}{A single character value specifying type of accession
IDs (\code{accession}). Must be specified when \code{accession} is specified.
(By default: \code{accession.type = NULL})}

\item{accession}{A single character value or a vector of character values
specifying accession IDs to return results for.
(By default: \code{accession = NULL})}

\item{as.df}{A single boolean value specifying whether to return the
results as a data.frame or leave as a nested list.
(By default: \code{as.df = TRUE})}
}
\value{
\code{data.frame} or \code{list}
}
\description{
Versatile function to retrieve raw results
}
\details{
This function returns data from MGnify database. Compared to
\code{getResult}, this function allows more flexible framework for fetching
the data. However, there are drawbacks: for counts data, \code{getResult}
returns optimally structured data container which is easier for downstream
analysis. \code{getData} returns raw data from the database. However, if
you want to retrieve data on pipelines or publications, for instance,
\code{getResult} is not suitable for it, and \code{getData} can be utilized
instead.
}
\examples{
# Create a client object
mg <- MgnifyClient(useCache = FALSE)

# Find kegg modules for certain analysis
df <- getData(
    mg, type = "kegg-modules",
    accession = "MGYA00642773", accession.type = "analyses")

}
\seealso{
\code{\link[MGnifyR:getResult]{getResult}}
}
