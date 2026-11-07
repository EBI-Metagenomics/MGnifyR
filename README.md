# MGnifyR <img src="man/figures//mgnifyr_logo.png" align="right" width="120" />

An R package for searching and retrieving data from the
[EBI Metagenomics resource](https://www.ebi.ac.uk/metagenomics) (MGnify).
Browse studies, and fetch their taxonomic profiles as a single
`TreeSummarizedExperiment`. The data is retrieved from the summary tables that
MGnify provides for each study, so fetching even hundreds of studies is fast.

MGnifyR package is part of [miaverse](https://microbiome.github.io/) 
microbiome analysis ecosystem enabling usage of
[mia](https://bioconductor.org/packages/release/bioc/html/mia.html)
and other miaverse packages.

<img src="man/figures/findingpheno_logo.png" align="right" width="160" />

**This research has received funding from the Horizon 2020 Programme of the
European Union within the FindingPheno project under grant agreement No
952914.** FindingPheno, an EU-funded project, is dedicated to developing
computational tools and methodologies for the integration and analysis of
multi-omics data. Its primary objective is to deepen our understanding of the
interactions between hosts and their microbiomes. You can find more information
on [FindingPheno website](https://findingpheno.eu/).

## Installation

### Bioc-release

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MGnifyR")
```

### Bioc-devel

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initialises usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("MGnifyR")
```

### GitHub

```
remotes::install_github("EBI-Metagenomics/MGnifyR")
```

## Basic usage

For more detailed instructions, read the function help and the vignette
(`vignette("MGnifyR")`).

```
library(MGnifyR)

# Browse studies
studies <- searchStudies(biome = "root:Host-associated:Human", max.hits = 10)

# Fetch their taxonomic profiles as a single TreeSummarizedExperiment
tse <- fetchStudies(studies)

# Or select specific analyses, e.g., based on sample metadata
analyses <- searchAnalyses("MGYS00005116")
tse <- fetchAnalyses(analyses[analyses$pipeline_version == "V5", ])
```
