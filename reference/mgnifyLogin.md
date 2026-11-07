# Log in to MGnify to access private data

Log in to MGnify to access private data

## Usage

``` r
mgnifyLogin(username, password)
```

## Arguments

- username:

  A single character value specifying the Webin username, e.g.,
  `"Webin-1122334"`.

- password:

  A single character value specifying the password.

## Value

`TRUE` invisibly, if login was successful.

## Details

Most MGnify data is public and can be fetched without logging in. To
access private data, log in with your ENA Webin credentials. The login
lasts for the R session, and it is used by all the functions that fetch
data from MGnify API v2, such as [`searchStudies`](searchStudies.md) and
[`fetchStudies`](fetchStudies.md).

## Examples

``` r
# Wrong credentials give an error
try(mgnifyLogin("Webin-0", "wrong"))
#> Error : Failed to log in.

# With your own ENA Webin credentials, log in and fetch private data
# mgnifyLogin("Webin-1122334", "SecretPassword")
# tse <- fetchStudies("MGYS00000001")
```
