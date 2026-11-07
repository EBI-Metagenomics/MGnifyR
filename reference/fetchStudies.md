# Fetch taxonomic and functional profiles

Fetch taxonomic and functional profiles

## Usage

``` r
fetchStudies(
  study.ids,
  type = "SSU",
  samples = FALSE,
  biosamples = FALSE,
  studies = FALSE,
  analyses = FALSE,
  assemblies = FALSE,
  qc = FALSE
)

fetchAnalyses(
  analysis.ids,
  type = "SSU",
  samples = FALSE,
  biosamples = FALSE,
  studies = FALSE,
  analyses = FALSE,
  assemblies = FALSE,
  qc = FALSE
)
```

## Arguments

- study.ids:

  A character vector specifying study accessions, e.g.,
  `"MGYS00005116"`, or a `data.frame` returned by
  [`searchStudies`](searchStudies.md).

- type:

  A character vector specifying the profiles to fetch, e.g., `"SSU"` or
  `c("SSU", "GO-slim")`. See details for the types. (By default:
  `type = "SSU"`)

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

- analysis.ids:

  A character vector specifying analysis accessions, e.g.,
  `"MGYA00383350"`, or a `data.frame` returned by
  [`searchAnalyses`](searchStudies.md). Giving the `data.frame` is
  faster.

## Value

For a single type, a `TreeSummarizedExperiment` with

- `assay`: `"counts"`

- `colData`: the columns of [`searchAnalyses`](searchStudies.md), i.e.,
  the metadata of the analyses, their studies and samples

- `rowData`: taxonomy ranks, or the descriptions of the rows of
  functional profiles

For multiple types, a `MultiAssayExperiment` with a
`TreeSummarizedExperiment` for each type that has data. The profiles are
linked by sample, e.g., the SSU profile of an amplicon analysis and the
GO profile of a metagenomic analysis of the same sample.

`NULL` is returned if no data could be fetched.

## Details

`fetchStudies` fetches all the analyses of the given studies, and
`fetchAnalyses` fetches specific analyses, e.g., a filtered result of
[`searchAnalyses`](searchStudies.md). Columns are analyses and rows are
taxa or, in functional profiles, functions, e.g., GO terms. Taxa are
matched across studies by their taxonomy lineage.

`type` tells which profile is fetched. Only analyses that have the type
are fetched.

- Taxonomic profiles: `"SSU"` and `"LSU"` (small and large subunit
  rRNA), `"mOTUs"` (marker genes of metagenomes and metatranscriptomes,
  pipelines v5-v6), `"UNITE"` and `"ITSoneDB"` (fungal ITS), `"PR2"`
  (eukaryotes), `"DADA2-SILVA"` and `"DADA2-PR2"` (amplicon sequence
  variants, pipeline v6), and `"contigs"` (number of contigs per taxon
  in assemblies, pipeline v6). Pipelines v1-v3 have a single profile,
  which is fetched with `"SSU"`.

- Functional profiles of metagenomes, metatranscriptomes and assemblies:
  `"InterPro"`, `"GO"` and `"GO-slim"` (pipelines v1-v6), `"Pfam"` and
  `"KO"` (KEGG orthologs, pipelines v5-v6), and `"antiSMASH"`
  (biosynthetic gene clusters of assemblies, pipelines v5-v6). In
  pipeline v6, only assemblies have functional profiles.

Profiles of different pipeline versions are based on different reference
databases, so they might not be directly comparable. To fetch a single
pipeline version, subset the result of
[`searchAnalyses`](searchStudies.md) by `pipeline_version`. Data is
cached for the R session; to keep it between sessions, set
`options(MGnifyR.cache = "path/to/cache")`.

The arguments `samples`, `biosamples`, `studies`, `analyses`,
`assemblies` and `qc` choose the metadata in `colData`, as in
[`searchAnalyses`](searchStudies.md). In `fetchAnalyses`, they are used
only if `analysis.ids` are accessions; for a `data.frame`, give them to
[`searchAnalyses`](searchStudies.md) instead.

## See also

[`searchStudies`](searchStudies.md)

## Examples

``` r
# Fetch all analyses of a study
tse <- fetchStudies("MGYS00005116")
#> Fetching 2 analyses one by one, because they are missing from the summary tables of studies. This can take a while.
#> The data includes analyses from multiple pipeline versions (V4.1, V5, V6). Their profiles are based on different reference databases. Consider subsetting the result of searchAnalyses() and fetching the analyses with fetchAnalyses().
tse
#> class: TreeSummarizedExperiment 
#> dim: 1014 6 
#> metadata(0):
#> assays(1): counts
#> rownames(1014): sk__Bacteria sk__Bacteria;k__;p__Acidobacteria ...
#>   sk__Eukaryota;k__Viridiplantae;p__Streptophyta;c__Magnoliopsida;o__Poales;f__Poaceae;g__Triticum
#>   sk__Eukaryota;k__Viridiplantae;p__Streptophyta;c__Magnoliopsida;o__Solanales;f__Solanaceae;g__Nicotiana;s__Nicotiana_attenuata
#> rowData names(8): superkingdom kingdom ... genus species
#> colnames(6): MGYA00383350 MGYA00383351 ... MGYA01026512 MGYA01026513
#> colData names(16): analysis_accession study_accession ...
#>   sample_updated_at abundance_table
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> rowLinks: NULL
#> rowTree: NULL
#> colLinks: NULL
#> colTree: NULL

# Select analyses and fetch them
analyses <- searchAnalyses("MGYS00005116")
analyses <- analyses[analyses[["pipeline_version"]] == "V5", ]
tse <- fetchAnalyses(analyses)

# \donttest{
# Fetch a functional profile of a study
tse <- fetchStudies("MGYS00005610", type = "GO-slim")
#> Fetching 2 analyses one by one, because they are missing from the summary tables of studies. This can take a while.

# Fetch taxonomic and functional profiles into MultiAssayExperiment
mae <- fetchStudies("MGYS00005610", type = c("SSU", "GO-slim"))
#> Fetching 2 analyses one by one, because they are missing from the summary tables of studies. This can take a while.

# Browse studies and fetch them
studies <- searchStudies(biome = "root:Host-associated:Human", max.hits = 5)
tse <- fetchStudies(studies)
#> Fetching 6 analyses one by one, because they are missing from the summary tables of studies. This can take a while.
#> Could not download 6 tables: 'https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_results/ERP001/ERP001506/version_1.0/ERS149129_A_SFF/otus/rdp_assigned_taxonomy/ERS149129_A_SFF_rRNAFiltered.fasta_rep_set_tax_assignments.txt', 'https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_results/ERP001/ERP001506/version_1.0/ERS149130_A_SFF/otus/rdp_assigned_taxonomy/ERS149130_A_SFF_rRNAFiltered.fasta_rep_set_tax_assignments.txt', 'https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_results/ERP001/ERP001506/version_1.0/ERS149131_A_SFF/otus/rdp_assigned_taxonomy/ERS149131_A_SFF_rRNAFiltered.fasta_rep_set_tax_assignments.txt', .... They might be missing from the server, or the server did not respond. If the server did not respond, try again later; only the missing tables are downloaded again.
#> No SSU abundance table could be fetched for 6 analyses (no such table, or it could not be downloaded): 'MGYA00002902', 'MGYA00002903', 'MGYA00002904', 'MGYA00002905', 'MGYA00002906', 'MGYA00002907'
#> The data includes analyses from multiple pipeline versions (V1, V2). Their profiles are based on different reference databases. Consider subsetting the result of searchAnalyses() and fetching the analyses with fetchAnalyses().
# }
```
