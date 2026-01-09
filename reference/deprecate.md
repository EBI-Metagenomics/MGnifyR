# These functions will be deprecated. Please use other functions instead.

These functions will be deprecated. Please use other functions instead.

## Usage

``` r
mgnify_client(
  username = NULL,
  password = NULL,
  usecache = FALSE,
  cache_dir = NULL,
  warnings = FALSE,
  use_memcache = FALSE,
  ...
)

mgnify_query(
  client,
  qtype = "samples",
  accession = NULL,
  asDataFrame = TRUE,
  maxhits = 200,
  usecache = FALSE,
  ...
)

mgnify_analyses_from_samples(client, accession, usecache = TRUE, ...)

mgnify_analyses_from_studies(client, accession, usecache = TRUE, ...)

mgnify_get_download_urls(
  client,
  accessions,
  accession_type,
  usecache = TRUE,
  ...
)

mgnify_download(
  client,
  url,
  file = NULL,
  read_func = NULL,
  usecache = TRUE,
  Debug = FALSE,
  ...
)

mgnify_get_analyses_results(
  client = NULL,
  accessions,
  retrievelist = c(),
  compact_results = TRUE,
  usecache = TRUE,
  bulk_dl = FALSE,
  ...
)

mgnify_get_analyses_phyloseq(
  client = NULL,
  accessions,
  usecache = TRUE,
  returnLists = FALSE,
  tax_SU = "SSU",
  get_tree = FALSE,
  ...
)

mgnify_get_analyses_metadata(client, accessions, usecache = TRUE, ...)

mgnify_retrieve_json(
  client,
  path = "biomes",
  complete_url = NULL,
  qopts = NULL,
  maxhits = 200,
  usecache = FALSE,
  Debug = FALSE,
  ...
)
```

## Arguments

- username:

  \-

- password:

  \-

- usecache:

  \-

- cache_dir:

  \-

- warnings:

  \-

- use_memcache:

  \-

- ...:

  \-

- client:

  \-

- qtype:

  \-

- accession:

  \-

- asDataFrame:

  \-

- maxhits:

  \-

- accessions:

  \-

- accession_type:

  \-

- url:

  \-

- file:

  \-

- read_func:

  \-

- Debug:

  \-

- retrievelist:

  \-

- compact_results:

  \-

- bulk_dl:

  \-

- returnLists:

  \-

- tax_SU:

  \-

- get_tree:

  \-

- path:

  \-

- complete_url:

  \-

- qopts:

  \-

## Value

\-
