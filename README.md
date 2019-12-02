# [MGnifyR](https://github.com/beadyallen/MGnifyR)

An R package for searching and retrieving data from the EBI Metagenomics resource. Currently undergoing heavy development, it nevertheless provides a useful set of tools to easily access and process MGnify data in R. In most cases, MGnifyR interacts directly with the JSONAPI, rather than relying on downloading analyses outputs as TSV files. Thus it is more general - allowing for example the intuitive combining of multiple studies and analyses into a single workflow, but is in some cases slower than the afformentioned direct access. Local caching of results on disk is implemented to help counter some of the overheads, but data downloads can be slow - particularly for functional annotation retrieval. 

## Requirements

```
devtools # for installation
phyloseq
plyr
dplyr
reshape2

httr
urltools
```

## Installation instructions
At the R terminal:
```
devtools::install_github("beadyallen/MGnifyR")
```


## Basic usage
For more detailed instructions read the associated function help and vignette (`vignette("MGNifyR")`)

```
library(MGnifyR)

#Set up the MGnify client instance
mgclnt <- mgnify_client(usecache = T, cache_dir = '/tmp/MGnify_cache')

#Retrieve the list of analyses associated with a study
accession_list <- mgnify_analyses_from_studies(mgclnt, "MGYS00005058", usecache = T)

#Download all associated study/sample and analysis metadata
meta_dataframe <- mgnify_get_analyses_metadata(mgclnt, accession_list, usecache = T )

#Convert analyses outputs to a single `phyloseq` object
psobj <- mgnify_get_analyses_phyloseq(mgclnt, meta_dataframe$analysis_accession, usecache = T)
psobj

#Retrieve Interpro assignment counts for these analyses
ip_df <- mgnify_get_analyses_results(mgclnt, meta_dataframe$analysis_accession, retrievelist = c("interpro-identifiers"), usecache = T)
head(ip_df)
```


