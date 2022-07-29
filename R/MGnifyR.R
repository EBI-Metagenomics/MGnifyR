# NOTES:

#
# Some useful keyboard shortcuts for package authoring:
#
#     Install Package: 'Ctrl + Shift + B'
#     Check Package: 'Ctrl + Shift + E'
#     Test Package: 'Ctrl + Shift + T'
#
# @import ape
# @import dplyr
# @import httr
# @import MultiAssayExperiment
# @import phyloseq
# @import plyr
# @import reshape2
# @import SingleCellExperiment
# @import SummarizedExperiment
# @import TreeSummarizedExperiment
# @import urltools


library(urltools)
library(phyloseq)
library(plyr)
library(reshape2)
library(dplyr)
library(ape)



##Example accessions:
# Analysis assembly: MGYA00379728
# Analysis metagenomic: MGYA00377678
# Analysis amplicon: MGYA00250889



baseurl <- 'https://www.ebi.ac.uk/metagenomics/api/v1'

##Filters possible - this comes from the django source code - would be nice if we could
# look it up.
# These DON'T seem to include all possible attributes ....
# And only some
sample_filters <- c('accession','experiment_type','biome_name','lineage','geo_loc_name','latitude_gte','latitude_lte',
                                     'longitude_gte','longitude_lte','species','instrument_model','instrument_platform','metadata_key',
                                     'metadata_value_gte','metadata_value_lte','metadata_value','environment_material','environment_feature',
                                     'study_accession','include')
biome_filters <- c('depth_gte','depth_lte')
study_filters <- c('accession','biome_name','lineage','centre_name','include')
run_filters <- c('accession','experiment_type','biome_name','lineage','species','instrument_platform','instrument_model',
                            # 'metadata_key','metadata_value_gte','metadata_value_lte','metadata_value','sample_accession','study_accession',
                            'include')
analysis_filters <- c('biome_name', 'lineage', 'experiment_type', 'species', 'sample_accession', 'pipeline_version')


#Combined together into a single queriably list
query_filters <- list(
    biomes <- biome_filters,
    samples <- sample_filters,
    studies <- study_filters,
    runs <- run_filters
)


### Result table caching

mgnify_memory_cache <- list()




#Internal functions to parse the attributes/hierarchy list into a data.frame
mgnify_parse_tax <- function(json){
    df <- as.data.frame(c(json$attributes["count"], unlist(json$attributes$hierarchy)), stringsAsFactors = F)
    df$index_id <- json$attributes$lineage
    df

}
mgnify_parse_func <- function(json){
    df <- as.data.frame(json$attributes, stringsAsFactors = F)
    df$index_id <- json$attributes$accession
    df
}

#Which parser do you use for which type of output?
#' @export
analyses_results_type_parsers <- list(taxonomy=mgnify_parse_tax,`taxonomy-itsonedb` = mgnify_parse_tax, `go-slim`=mgnify_parse_func,
                                                                            `taxonomy-itsunite`=mgnify_parse_tax, `taxonomy-ssu`=mgnify_parse_tax,
                                                                            `taxonomy-lsu`=mgnify_parse_tax,`antismash-gene-clusters`=mgnify_parse_func,
                                                                            `go-terms`=mgnify_parse_func, `interpro-identifiers`=mgnify_parse_func)


#this maps the json attribute name for retrievelist to the "description -> label" attribute in the study downloads section
analyses_results_bulk_file_names <- list(
                                                                        #taxonomy="Taxonomic assignments SSU",
                                                                        `taxonomy-itsonedb` = "Taxonomic assignments ITS",
                                                                        `go-slim`="GO slim annotation",
                                                                        `taxonomy-itsunite`="Taxonomic assignments Unite", `taxonomy-ssu`="Taxonomic assignments SSU",
                                                                        `taxonomy-lsu`="Taxonomic assignments LSU",`antismash-gene-clusters`=mgnify_parse_func,
                                                                        `go-terms`="Complete GO annotation", `interpro-identifiers`="InterPro matches",
                                                                        `phylo-tax-ssu`="Phylum level taxonomies SSU",`phylo-tax-lsu`="Phylum level taxonomies LSU" )



##' @exportClass mgnify_client
mgnify_client <- setClass("mgnify_client",slots=list(url = "character", authtok="character",cache_dir="character",warnings="logical",use_memcache="logical",memcache="list",clear_cache="logical"),prototype = list(url=baseurl, authtok=NULL, cache_dir=NULL, use_memcache=FALSE, memcache=list(),clear_cache=FALSE))

