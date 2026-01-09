# Look up analysis accession IDs for one or more study or sample accessions

Look up analysis accession IDs for one or more study or sample
accessions

## Usage

``` r
searchAnalysis(x, ...)

# S4 method for class 'MgnifyClient'
searchAnalysis(x, type, accession, ...)
```

## Arguments

- x:

  A `MgnifyClient` object.

- ...:

  Optional arguments; not currently used.

- type:

  A single character value specifying a type of accession IDs specified
  by `accession`. Must be "studies" or "samples".

- accession:

  A single character value or a vector of character values specifying
  study or sample accession IDs that are used to retrieve analyses IDs.

## Value

Vector of analysis accession IDs.

## Details

Retrieve analysis accession IDs associated with the supplied study or
sample accession. In MGnify, an analysis accession refers to a certain
pipeline analysis, such as specific 16S rRNA or shotgun metagenomic
mapping. Studies can include multiple samples, and each sample can
undergo multiple analyses using these pipelines. Each analysis is
identified by a unique accession ID, allowing precise tracking and
retrieval of analysis results within the MGnify database.

## Examples

``` r
# Create a client object
mg <- MgnifyClient(useCache = FALSE)

# Retrieve analysis ids from study MGYS00005058
result <- searchAnalysis(mg, "studies", c("MGYS00005058"))
#> Fetching analyses...
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%

if (FALSE) { # \dontrun{
# Retrieve all analysis ids from samples
result <- searchAnalysis(
    mg, "samples", c("SRS4392730", "SRS4392743"))
} # }
```
