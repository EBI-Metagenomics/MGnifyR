# Constructor for creating a MgnifyClient object to allow the access to MGnify database.

Constructor for creating a MgnifyClient object to allow the access to
MGnify database.

A MgnifyClient object

## Usage

``` r
MgnifyClient(
  username = NULL,
  password = NULL,
  useCache = FALSE,
  cacheDir = tempdir(),
  showWarnings = FALSE,
  verbose = TRUE,
  clearCache = FALSE,
  ...
)
```

## Arguments

- username:

  A single character value specifying an optional username for
  authentication. (By default: `username = NULL`)

- password:

  A single character value specifying an optional password for
  authentication. (By default: `password = NULL`)

- useCache:

  A single boolean value specifying whether to enable on-disk caching of
  results during this session. In most use cases should be TRUE. (By
  default: `useCache = FALSE`)

- cacheDir:

  A single character value specifying a folder to contain the local
  cache. Note that cached files are persistent, so the cache directory
  may be reused between sessions, taking advantage of previously
  downloaded results. The directory will be created if it doesn't exist
  already. (By default: `cacheDir = tempdir()`)

- showWarnings:

  A single boolean value specifying whether to print warnings during
  invocation of some MGnifyR functions. (By default:
  `showWarnings = FALSE`)

- verbose:

  A single boolean value specifying whether to print extra output during
  invocation of some MGnifyR functions. (By default: `verbose = FALSE`)

- clearCache:

  A single boolean value specifying whether to clear the cache. (By
  default: `clearCache = FALSE`)

- ...:

  optional arguments:

  - **url** A single character value specifying an url address of the
    database. (By default:
    `url = "https://www.ebi.ac.uk/metagenomics/api/v1"`)

## Value

A MgnifyClient object.

## Details

All functions in the MGnifyR package take a `MgnifyClient` object as
their first argument. The object allows the simple handling of both user
authentication and access to private data, and manages general options
for querying the MGnify database.

An object that are required by functions of MGnifyR package.

## Slots

- `databaseUrl`:

  A single character value specifying an URL address of database.

- `authTok`:

  A single character value specifying authentication token.

- `useCache`:

  A single boolean value specifying whether to use cache.

- `cacheDir`:

  A single character value specifying cache directory.

- `showWarnings`:

  A single boolean value specifying whether to show warnings.

- `clearCache`:

  A single boolean value specifying whether to clear cache.

- `verbose`:

  A single boolean value specifying whether to show messages.

## Constructor

See `MgnifyClient` for constructor.

## Accessor

See [`MgnifyClient-accessors`](MgnifyClient-accessors.md) for accessor
functions.

## Examples

``` r
my_client <- MgnifyClient(
    useCache = TRUE, cacheDir = "/scratch/MGnify_cache_location"
    )

if (FALSE) { # \dontrun{
# Use username and password to get access to non-public data
my_client <- MgnifyClient(
    username = "Webin-1122334", password = "SecretPassword",
    useCache = TRUE, cacheDir = "/scratch/MGnify_cache_location"
    )
} # }
```
