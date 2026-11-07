# Deprecated and defunct functions

Deprecated and defunct functions

## Usage

``` r
MgnifyClient(
  username = NULL,
  password = NULL,
  useCache = FALSE,
  cacheDir = tempdir(),
  ...
)

getResult(
  x,
  accession,
  get.taxa = "SSU",
  get.func = TRUE,
  output = "TreeSE",
  ...
)

searchAnalysis(x, type, accession, ...)

getMetadata(x, accession, ...)

doQuery(
  x,
  type = "studies",
  accession = NULL,
  as.df = TRUE,
  max.hits = 200,
  ...
)

getData(x, type, accession.type = NULL, accession = NULL, as.df = TRUE, ...)

searchFile(
  x,
  accession,
  type = c("studies", "samples", "analyses", "assemblies", "genomes", "run"),
  ...
)

getFile(x, url, file = NULL, read.func = NULL, ...)

databaseUrl(...)

authTok(...)

useCache(...)

cacheDir(...)

showWarnings(...)

clearCache(...)

verbose(...)

databaseUrl(x) <- value

authTok(x) <- value

useCache(x) <- value

cacheDir(x) <- value

showWarnings(x) <- value

clearCache(x) <- value

verbose(x) <- value

mgnify_client(...)

mgnify_query(...)

mgnify_analyses_from_samples(...)

mgnify_analyses_from_studies(...)

mgnify_get_download_urls(...)

mgnify_download(...)

mgnify_get_analyses_results(...)

mgnify_get_analyses_phyloseq(...)

mgnify_get_analyses_metadata(...)

mgnify_retrieve_json(...)
```

## Arguments

- username, password:

  Webin credentials for [`mgnifyLogin`](mgnifyLogin.md).

- useCache, cacheDir:

  If `useCache = TRUE`, `cacheDir` is set as the cache directory.

- ...:

  In `doQuery`, query parameters of the API. Otherwise ignored.

- x:

  Ignored. Formerly a `MgnifyClient` object.

- accession:

  A character vector of accessions.

- get.taxa:

  Taxonomy to fetch, e.g., `"SSU"`.

- get.func, output:

  Not supported anymore.

- type:

  The type of accessions in `searchAnalysis` and `searchFile`, or the
  API endpoint in `doQuery` and `getData`, e.g., `"studies"`.

- as.df:

  Whether to return a `data.frame` instead of a list.

- max.hits:

  The maximum number of results, or `NULL` for all.

- accession.type:

  The type of accessions in `getData`.

- url:

  The URL of the file.

- file:

  `NULL` or a path to save the file to.

- read.func:

  `NULL` or a function applied to the downloaded file.

- value:

  Ignored.

## Value

See the replacing functions. `doQuery`, `getData` and `searchFile`
return a `data.frame`, and `getFile` returns the path of the file or the
result of `read.func`.

## Details

MGnify API v1 is frozen, and MGnifyR now uses MGnify API v2. The
functions of the earlier versions of MGnifyR are deprecated, and they
will be removed in a future release. Where possible, they use the new
functions:

- `MgnifyClient`: a client is not needed anymore. Log in with
  [`mgnifyLogin`](mgnifyLogin.md) to access private data, and set
  `options(MGnifyR.cache = "path")` to keep the cache between sessions.
  For compatibility, `MgnifyClient` does these, and returns `NULL`. The
  `x` argument of the other deprecated functions is ignored.

- `getResult`: use [`fetchAnalyses`](fetchStudies.md). Only taxonomic
  profiles are returned, as `TreeSummarizedExperiment`.

- `searchAnalysis` and `getMetadata`: use
  [`searchAnalyses`](searchStudies.md).

- `doQuery` and `getData`: use [`searchStudies`](searchStudies.md),
  [`searchSamples`](searchStudies.md) or
  [`searchAnalyses`](searchStudies.md). These query any endpoint of
  [MGnify API v2](https://www.ebi.ac.uk/metagenomics/api/v2/), but the
  endpoints and filters differ from API v1.

- `searchFile` and `getFile`: list and download the files of studies,
  samples or analyses.

The accessors of `MgnifyClient`, such as `useCache`, and the functions
with `mgnify_` prefix are defunct.

## Examples

``` r
# Instead of getResult(), use fetchAnalyses()
tse <- fetchAnalyses("MGYA00383350")
```
