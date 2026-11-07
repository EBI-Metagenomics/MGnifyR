# Search MGnify studies, samples and analyses

Search MGnify studies, samples and analyses

## Usage

``` r
searchStudies(
  biome = NULL,
  pipeline.version = NULL,
  search = NULL,
  order = NULL,
  max.hits = NULL
)

searchSamples(
  study.ids = NULL,
  biome = NULL,
  search = NULL,
  order = NULL,
  max.hits = NULL
)

searchAnalyses(
  study.ids,
  samples = FALSE,
  biosamples = FALSE,
  studies = FALSE,
  analyses = FALSE,
  assemblies = FALSE,
  qc = FALSE
)
```

## Arguments

- biome:

  `NULL` or a single character value specifying a biome lineage, e.g.,
  `"root:Environmental:Aquatic"`. Results also include the descendant
  biomes. Matching is not case-sensitive. The available biomes are
  listed at <https://www.ebi.ac.uk/metagenomics/browse/biomes>. (By
  default: `biome = NULL`)

- pipeline.version:

  `NULL` or a single character value specifying the version of the
  MGnify pipeline: `"V1"`, `"V2"`, `"V3"`, `"V4"`, `"V4.1"`, `"V5"`,
  `"V6"`, `"V6.1"` or `"V6.2"`. The `"V"` and `".0"` can be omitted,
  e.g., `"5.0"` or `"5"`. (By default: `pipeline.version = NULL`)

- search:

  `NULL` or a single character value used to search titles and
  accessions. (By default: `search = NULL`)

- order:

  `NULL` or a single character value specifying the order of the
  results. For studies, `"accession"` or `"updated_at"`, and for
  samples, `"sample_title"` or `"updated_at"`. Prefix with `"-"` for
  descending order, e.g., `"-updated_at"`. (By default: `order = NULL`)

- max.hits:

  `NULL` or a single integer value specifying the maximum number of
  results. If `NULL`, all results are returned. (By default:
  `max.hits = NULL`)

- study.ids:

  A character vector specifying study accessions, e.g.,
  `"MGYS00005116"`, or a `data.frame` returned by `searchStudies`. In
  `searchSamples`, it can also be `NULL` to search all samples. `biome`,
  `search` and `order` can be used only when searching all samples.

- samples:

  `TRUE` or `FALSE` specifying whether to add the sample metadata in
  MGnify. (By default: `samples = FALSE`)

- biosamples:

  `TRUE` or `FALSE` specifying whether to add the attributes of the
  samples in BioSamples. (By default: `biosamples = FALSE`)

- studies:

  `TRUE` or `FALSE` specifying whether to add the details and
  publications of the studies. (By default: `studies = FALSE`)

- analyses:

  `TRUE` or `FALSE` specifying whether to add the details of the
  analyses, e.g., the quality control of pipeline v6. (By default:
  `analyses = FALSE`)

- assemblies:

  `TRUE` or `FALSE` specifying whether to add the details of the
  assemblies. (By default: `assemblies = FALSE`)

- qc:

  `TRUE` or `FALSE` specifying whether to add the quality control of
  pipelines v1-v5, e.g., the number of sequences before and after
  quality filtering. It is read from separate files, which makes the
  search slower. The details of the analyses are added too. (By default:
  `qc = FALSE`)

## Value

A `data.frame` where each row is a study, sample or analysis.

## Details

Browse studies with `searchStudies` and fetch their data with
[`fetchStudies`](fetchStudies.md). To select specific analyses, e.g., by
sample metadata or pipeline version, filter the result of
`searchAnalyses` and give it to [`fetchAnalyses`](fetchStudies.md).

`searchAnalyses` returns one row per analysis with the accessions,
pipeline version, experiment type, sequencing instrument, and the other
fields of the analysis in the list of analyses of its study. The other
sources are added with separate arguments, as each needs its own
requests. By default, none of them is added, which is the fastest. The
data is cached, so a later search with more sources fetches only the new
data:

- `samples`: sample metadata in MGnify (`sample_*`), one request per 100
  samples

- `biosamples`: all attributes of the samples in
  [BioSamples](https://www.ebi.ac.uk/biosamples/) (`biosamples_*`),
  including study-specific variables, e.g., `biosamples_disease_status`,
  one request per 100 samples

- `studies`: details and publications of the studies (`study_*`), two
  requests per study

- `analyses`: details of the analyses (`analysis_*`), e.g., the quality
  control of pipeline v6, one request per analysis

- `assemblies`: details of the assemblies (`assembly_*`), one request
  per assembly

- `qc`: quality control files of pipelines v1-v5 (`analysis_*`), several
  files per analysis, and the details of the analyses

These are the same columns as in the `colData` of
[`fetchAnalyses`](fetchStudies.md). Values are as in the database.
Lists, e.g., the publications of a study (`study_publications`), are
JSON text that can be read with
[`jsonlite::fromJSON()`](https://jeroen.r-universe.dev/jsonlite/reference/fromJSON.html).

In MGnify, a study has samples, a sample has runs, i.e., sequencings of
the sample, and an analysis is a run or an assembly processed with an
MGnify pipeline. A sample can have multiple runs, e.g., an amplicon and
a metagenome, and a run can be analysed with multiple pipeline versions,
so a sample can have multiple analyses. An analysis usually has one run,
but some have multiple runs, e.g., `"ERR1, ERR2"`, and analyses of
assemblies can have no run in `run_accession`.

Results are cached for the R session; to keep them between sessions, set
`options(MGnifyR.cache = "path/to/cache")`. To access private data, log
in first with [`mgnifyLogin`](mgnifyLogin.md).

## See also

[`fetchStudies`](fetchStudies.md)

## Examples

``` r
# Search studies of aquatic environments
studies <- searchStudies(biome = "root:Environmental:Aquatic", max.hits = 10)

# Get samples and analyses of a study
samples <- searchSamples("MGYS00005116")
analyses <- searchAnalyses("MGYS00005116")
analyses <- analyses[analyses[["pipeline_version"]] == "V6", ]
```
