# MGnifyR <img src="man/figures//mgnifyr_logo.png" align="right" width="120" />

An R package for searching and retrieving data from the EBI Metagenomics resource. 
In most cases, MGnifyR interacts directly with the JSONAPI, rather than relying
on downloading analyses outputs as TSV files. Thus it is more general - allowing
for example the intuitive combining of multiple studies and analyses
into a single workflow, but is in some cases slower than the afformentioned
direct access. Local caching of results on disk is implemented to help counter
some of the overheads, but data downloads can be slow - particularly for
functional annotation retrieval. 

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

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("MGnifyR")
```

### GitHub

```
remotes::install_github("EBI-Metagenomics/MGnifyR")
```

## Basic usage
For more detailed instructions read the associated function help and vignette (`vignette("MGNifyR")`)

```
library(MGnifyR)

# Set up the MGnify client instance
mgclnt <- MgnifyClient(usecache = TRUE, cache_dir = '/tmp/MGnify_cache')

# Retrieve the list of analyses associated with a study
accession_list <- searchAnalysis(mgclnt, "studies", "MGYS00005058", usecache = TRUE)

# Download all associated study/sample and analysis metadata
meta_dataframe <- getMetadata(mgclnt, accession_list, usecache = TRUE)

# Convert analyses outputs to a single `MultiAssayExperiment` object
mae <- getResult(mgclnt, meta_dataframe$analysis_accession, usecache = TRUE)
mae
```

