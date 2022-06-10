devtools::load_all()

library(MGnifyR)

# Set up the MGnify client instance

mgclnt <- mgnify_client(usecache = T, cache_dir = '/tmp/MGnify_cache')

# Retrieve the list of analyses associated with a study

accession_list <- mgnify_analyses_from_studies(mgclnt, "MGYS00005058", usecache = T)

# Download all associated study/sample and analysis metadata

meta_dataframe <- mgnify_get_analyses_metadata(mgclnt, accession_list, usecache = T )

# Convert analyses outputs to a single `phyloseq` object

psobj <- mgnify_get_analyses_phyloseq(mgclnt, meta_dataframe$analysis_accession, usecache = T)

######################## kev part #####################################
tse <- mgnify_get_analyses_tse(mgclnt, meta_dataframe$analysis_accession, usecache = T)
tse
###################### end kev part ###################################

# retrievelist requires us to specify which functional estimations we

# wish to examine, and may be one or more of go-slim, go-terms,

# interpro-identifiers or antismash-gene-clusters. retrievelist may

# also be all, in which case all available results are retrieved,

# along with all taxonomic assignments.

# Retrieve assignment counts for these analyses as a data list

#ip_df <- mgnify_get_analyses_results(mgclnt, meta_dataframe$analysis_accession, retrievelist = c("interpro-identifiers"), usecache = T)

dli <- mgnify_get_analyses_results(mgclnt, meta_dataframe$analysis_accession, retrievelist = "all", usecache = T)
dim(dli)
length(dli)

# Seems that taxonomy-ssu corresponds to the phyloseq object

# all(abundances(psobj) == dli[["taxonomy-ssu"]][, 9:32])

# TODO

# 1. mgnify_get_analyses_phyloseq -> mgnify_get_analyses_treese (new function)

# 2. mgnify_get_analyses_results -> mgnify_get_analyses_results_mae (new function; to download & arrange selected tables in a MultiAssayExperiment)

# 3. some use examples in MGnifyR vignette

# 4. some use examples in OMA

# 5. MGnifyR to CRAN (joint project?)

