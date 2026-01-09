# Download any MGnify files, also including processed reads and identified protein sequences

Download any MGnify files, also including processed reads and identified
protein sequences

Listing files available for download

## Usage

``` r
getFile(x, ...)

searchFile(x, ...)

# S4 method for class 'MgnifyClient'
getFile(x, url, file = NULL, read.func = NULL, ...)

# S4 method for class 'MgnifyClient'
searchFile(
  x,
  accession,
  type = c("studies", "samples", "analyses", "assemblies", "genomes", "run"),
  ...
)
```

## Arguments

- x:

  A `MgnifyClient` object.

- ...:

  Additional arguments; not used currently.

- url:

  A single character value specifying the url address of the file we
  wish to download.

- file:

  A single character value or NULL specifying an optional local filename
  to use for saving the file. If `NULL`, MGNify local cache settings
  will be used. If the file is intended to be processed in a separate
  program, it may be sensible to provide a meaningful `file`, rather
  than having to hunt through the cache folders. If `file` is `NULL` and
  `useCache(client)` is `FALSE`, the `read.func` parameter must be
  supplied or the file will be downloaded and then deleted. (By default:
  `file = NULL`)

- read.func:

  A function specifying an optional function to process the downloaded
  file and return the results, rather than relying on post processing.
  The primary use-case for this parameter is when local disk space is
  limited and downloaded files can be quickly processed and discarded.
  The function should take a single parameter, the downloaded filename,
  and may return any valid R object. (By default: `read.func = NULL`)

- accession:

  A single character value or a vector of character values specifying
  accession IDs to return results for.

- type:

  A single character value specifying the type of objects to query. Must
  be one of the following options: `analysis`, `samples`, `studies`,
  `assembly`, `genome` or `run`. (By default: `type = "samples"`)

## Value

For `getFile()`, either the local filename of the downloaded file, be it
either the location in the MGNifyR cache or file. If `read.func` is
used, its result will be returned.

For `searchFile()` `data.frame` containing all discovered downloads. If
multiple `accessions` are queried, the `accessions` column may to filter
the results - since rownames are not set (and wouldn't make sense as
each query will return multiple items)

## Details

`getFile` is a convenient wrapper round generic the URL downloading
functionality in R, taking care of things like local caching and
authentication.

`searchFile()` function is a wrapper function allowing easy enumeration
of downloads available for a given accession IDs. Returns a single
data.frame containing all available downloads and associated metadata,
including the url location and description. This can then be filtered to
extract the urls of interest, before actually retrieving the files using
`getFile()`

## Examples

``` r
# Make a client object
mg <- MgnifyClient(useCache = FALSE)

# Create a vector of accession ids - these happen to be \code{analysis}
# accessions
accession_vect <- c("MGYA00563876", "MGYA00563877")
downloads <- searchFile(mg, accession_vect, "analyses")
#> Searching files...
#>   |                                                                              |                                                                      |   0%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%

# Filter to find the urls of 16S encoding sequences
url_list <- downloads[
    downloads$attributes.description.label == "Contigs encoding SSU rRNA",
    "download_url"]

# Example 1:
# Download the first file
supplied_filename <- getFile(
    mg, url_list[[1]], file="SSU_file.fasta.gz")

if (FALSE) { # \dontrun{
# Example 2:
# Just use local caching
cached_filename <- getFile(mg, url_list[[2]])

# Example 3:
# Using read.func to open the reads with readDNAStringSet from
# \code{biostrings}. Without retaining on disk
dna_seqs <- getFile(
    mg, url_list[[3]], read.func = readDNAStringSet)
} # }

# Make a client object
mg <- MgnifyClient(useCache = TRUE)
# Create a vector of accession ids - these happen to be \code{analysis}
# accessions
accession_vect <- c(
    "MGYA00563876", "MGYA00563877", "MGYA00563878",
    "MGYA00563879", "MGYA00563880" )
downloads <- searchFile(mg, accession_vect, "analyses")
#> Searching files...
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
```
