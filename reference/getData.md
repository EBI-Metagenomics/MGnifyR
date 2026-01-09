# Versatile function to retrieve raw results

Versatile function to retrieve raw results

## Usage

``` r
getData(x, ...)

# S4 method for class 'MgnifyClient'
getData(x, type, accession.type = NULL, accession = NULL, as.df = TRUE, ...)
```

## Arguments

- x:

  A `MgnifyClient` object.

- ...:

  optional arguments fed to internal functions.

- type:

  A single character value specifying the type of data retrieve. Must be
  one of the following options: `studies`, `samples`, `runs`,
  `analyses`, `biomes`, `assemblies`, `super-studies`,
  `experiment-types`, `pipelines`, `pipeline-tools`, `publications`,
  `genomes`, `genome-search`, `genome-search/gather`,
  `genome-catalogues`, `genomeset`, `cogs`, `kegg-modules`,
  `kegg-classes`, `antismash-geneclusters`, `annotations/go-terms`,
  `annotations/interpro-identifiers`, `annotations/kegg-modules`,
  `annotations/pfam-entries`, `annotations/kegg-orthologs`,
  `annotations/genome-properties`,
  `annotations/antismash-gene-clusters`, `annotations/organisms`, or
  `mydata`.

- accession.type:

  A single character value specifying type of accession IDs
  (`accession`). Must be specified when `accession` is specified. (By
  default: `accession.type = NULL`)

- accession:

  A single character value or a vector of character values specifying
  accession IDs to return results for. (By default: `accession = NULL`)

- as.df:

  A single boolean value specifying whether to return the results as a
  data.frame or leave as a nested list. (By default: `as.df = TRUE`)

## Value

`data.frame` or `list`

## Details

This function returns data from MGnify database. Compared to
`getResult`, this function allows more flexible framework for fetching
the data. However, there are drawbacks: for counts data, `getResult`
returns optimally structured data container which is easier for
downstream analysis. `getData` returns raw data from the database.
However, if you want to retrieve data on pipelines or publications, for
instance, `getResult` is not suitable for it, and `getData` can be
utilized instead.

## See also

[`getResult`](getResult.md)

## Examples

``` r
# Create a client object
mg <- MgnifyClient(useCache = FALSE)

# Find kegg modules for certain analysis
df <- getData(
    mg, type = "kegg-modules",
    accession = "MGYA00642773", accession.type = "analyses")
```
