\name{mgnify_query_samples}
\alias{mgnify_query_samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
test123
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mgnify_query_samples(client, accession = NULL, asDataFrame = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{client}{
%%     ~~Describe \code{client} here~~
}
  \item{accession}{
%%     ~~Describe \code{accession} here~~
}
  \item{asDataFrame}{
%%     ~~Describe \code{asDataFrame} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (client, accession = NULL, asDataFrame = F, ...)
{
    a = accession
    arglist = as.list(match.call())
    arglist$accession = a
    qopts = arglist[names(arglist) \%in\% sample_filters]
    result = mgnify_query_json(client, path = "samples", qopts = qopts)
    samp_id_list = lapply(result, function(x) x$id)
    names(result) = samp_id_list
    if (asDataFrame) {
        samplist = sapply(result, function(r) {
            attrlist = names(r$attributes)
            baseattrlist = attrlist[!(attrlist \%in\% "sample-metadata")]
            metaattrlist = r$attributes$`sample-metadata`
            metlist = sapply(metaattrlist, function(x) x$value)
            names(metlist) = sapply(metaattrlist, function(x) x$key)
            df = as.data.frame(t(unlist(c(r$attributes[baseattrlist],
                metlist))), stringsAsFactors = F)
            df$biome = r$relationship$biome$data$id
            df$study = r$relationship$studies$data[[1]]$id
            rownames(df) = df$accession
            df
        })
        tryCatch(plyr::rbind.fill(samplist), error = function(e) samplist)
    }
    else {
        result
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
