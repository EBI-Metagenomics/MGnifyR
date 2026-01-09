# Get all study, sample and analysis metadata for the supplied analysis accessions

Get all study, sample and analysis metadata for the supplied analysis
accessions

## Usage

``` r
getMetadata(x, ...)

# S4 method for class 'MgnifyClient'
getMetadata(x, accession, ...)
```

## Arguments

- x:

  A `MgnifyClient` object.

- ...:

  Optional arguments; not currently used.

- accession:

  A single character value or a vector of analysis accession IDs
  specifying accessions to retrieve data for.

## Value

A `data.frame` containing metadata for each analysis in the `accession`
list. Each row represents a single analysis.

## Details

The function retrieves all study, sample and analysis metadata
associated with provided analysis accessions.

## Examples

``` r
# Create a client object
mg <- MgnifyClient(useCache = FALSE)

# Download all associated study/sample and analysis metadata
accession_list <- c("MGYA00377505")
meta_dataframe <- getMetadata(mg, accession_list)
#> Fetching metadata...
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
```
