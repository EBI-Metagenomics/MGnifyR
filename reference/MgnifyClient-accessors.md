# MgnifyClient accessors and mutators

MgnifyClient accessors and mutators

## Usage

``` r
databaseUrl(x)

authTok(x)

useCache(x)

cacheDir(x)

showWarnings(x)

clearCache(x)

verbose(x)

databaseUrl(x) <- value

authTok(x) <- value

useCache(x) <- value

cacheDir(x) <- value

showWarnings(x) <- value

clearCache(x) <- value

verbose(x) <- value

# S4 method for class 'MgnifyClient'
databaseUrl(x)

# S4 method for class 'MgnifyClient'
authTok(x)

# S4 method for class 'MgnifyClient'
useCache(x)

# S4 method for class 'MgnifyClient'
cacheDir(x)

# S4 method for class 'MgnifyClient'
showWarnings(x)

# S4 method for class 'MgnifyClient'
clearCache(x)

# S4 method for class 'MgnifyClient'
verbose(x)

# S4 method for class 'MgnifyClient'
databaseUrl(x) <- value

# S4 method for class 'MgnifyClient'
authTok(x) <- value

# S4 method for class 'MgnifyClient'
useCache(x) <- value

# S4 method for class 'MgnifyClient'
cacheDir(x) <- value

# S4 method for class 'MgnifyClient'
showWarnings(x) <- value

# S4 method for class 'MgnifyClient'
clearCache(x) <- value

# S4 method for class 'MgnifyClient'
verbose(x) <- value
```

## Arguments

- x:

  A `MgnifyClient` object.

- value:

  A value to be added to a certain slot.

## Value

A value of MgnifyClient object or nothing.

## Details

These functions are for fetching and mutating slots of `MgnifyClient`
object.

## Examples

``` r
mg <- MgnifyClient()

databaseUrl(mg)
#> [1] "https://www.ebi.ac.uk/metagenomics/api/v1"
showWarnings(mg) <- FALSE
```
