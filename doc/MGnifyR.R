## ---- include = FALSE----------------------------------------------------
library(knitr)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache = TRUE
)

## ----devtools_install, eval=FALSE----------------------------------------
#  devtools::install_github("beadyallen/MGnifyR")

## ----load_package--------------------------------------------------------
library(MGnifyR)

## ----create_client, echo=TRUE, fig.keep='all', message = FALSE-----------
mg <- mgnify_client( usecache=T, cache_dir='/home/ben/.MGnify_cache')

## ----create_client_passwd, eval=FALSE------------------------------------
#  mg <- mgnify_client(username="Webin-username", password="your-password", usecache = T)

## ----search_studies, fig.keep='all', message = FALSE---------------------
northpolar <- mgnify_query(mg, "samples", latitude_gte=60.0, experiment_type="amplicon", biome_name="Soil", instrument_platform = "Illumina", usecache = F )
northpolar[1:5,]

## ----search_studies_accession, message = FALSE---------------------------
study_samples <- mgnify_query(mg, "samples", study_accession="MGYS00003725", usecache=T)
#kable(study_samples[1:10,])
study_samples[1:5,]

## ----convert_to_analyses, fig.keep='all', results='hide', message = FALSE----
#just retrieve 20 hits for demonstration
analyses_accessions <- mgnify_analyses_from_samples(mg, accession = study_samples$accession[1:20])

## ----show_accessions-----------------------------------------------------
analyses_accessions

## ----get_metadata, fig.keep='all', results='hide', message = FALSE-------
metadata <- mgnify_get_analyses_metadata(mg, analyses_accessions)

## ----show_metadata-------------------------------------------------------
head(metadata)

## ----get_analyses, results='hide', fig.keep='all', message = FALSE-------

soil <- mgnify_analyses_from_studies(mg, "MGYS00001447")[1:20]
human <- mgnify_analyses_from_studies(mg, "MGYS00001442")[1:20]
marine <- mgnify_analyses_from_studies(mg, "MGYS00001282")[1:20]

all_accessions <- c(soil,human,marine)


## ----get__new_metadata, echo=TRUE, results='hide', fig.keep='all', message = FALSE----
full_metadata <- mgnify_get_analyses_metadata(mg, all_accessions)

## ----show_new_metadata---------------------------------------------------
colnames(full_metadata)
head(full_metadata)

## ----full_metatdata_explore----------------------------------------------
#Distribution of sample source material:
table(full_metadata$`sample_environment-material`)

#What sequencing machine(s) were used?
table(full_metadata$`sample_instrument model`)

# Boxplot of raw read counts:
ggplot(full_metadata, aes(x=study_accession, y=log(as.numeric(`analysis_Submitted nucleotide sequences`)))) + geom_boxplot(aes(group=study_accession)) + theme_bw() + ylab("log(submitted reads)")


## ----get_phyloseq, echo=TRUE, results='hide', fig.keep='all', message = FALSE----
full_phyloseq <- mgnify_get_analyses_phyloseq(mg, full_metadata$analysis_accession)

## ----show_phyloseq-------------------------------------------------------
full_phyloseq

## ----plot_taxa, echo=TRUE, fig.align="center", fig.height=4, fig.width=6, fig.keep='all', message = FALSE----
library(phyloseq)
library(ggplot2)

#rarefy the data (~sorry~) for alpha diversity
normed_ps <- rarefy_even_depth(full_phyloseq, rngseed=1)

class_ps <- tax_glom(normed_ps, "Class")
plot_bar(class_ps,  fill="Phylum")  + theme_bw() + theme(legend.position = "none")

## ----plot_diversity, fig.keep='all', message = FALSE, fig.align="center", fig.height=4, fig.width=6----
alphadiversity = estimate_richness(normed_ps)

adf <- cbind.data.frame(phyloseq::sample_data(normed_ps)$`sample_environment.biome`, alphadiversity$InvSimpson)
colnames(adf) <- c("study","diversity")
ggplot(adf, aes(x=study, group=study, y=diversity)) + geom_boxplot() + theme_bw()

## ---- get_functions, cache=TRUE, echo=TRUE, results='hide', fig.keep='all', message = FALSE----

func_res <- mgnify_get_analyses_results(mg, full_metadata$analysis_accession, retrievelist = "go-slim", bulk_dl = T)

goslim <- func_res$`go-slim`

## ----show_goslim---------------------------------------------------------
head(goslim)


## ---- nmds_function, results='hide', fig.keep='all', message = FALSE, fig.align="center", fig.height=4, fig.width=6----
library(vegan)

#Find the per-sample raw read count for each sample
seqvect <- as.numeric(full_metadata[colnames(goslim)[-c(1,2,3)],"analysis_Nucleotide sequences after format-specific filtering"])

#scale factors for normalizing the go term results
scale_factors = 1/(seqvect/median(seqvect))

#Extract the numeric matrix
normalized= goslim[,-c(1,2,3)] * scale_factors

#Calculate NMDS 
nmds_res = vegan::metaMDS(t(normalized))

#Build the plot
results_df <-  merge(nmds_res$points, 
                     full_metadata[,c("analysis_accession","sample_environment-feature")], 
                     by.x="row.names", by.y="analysis_accession")

ggplot(results_df, aes(x=MDS1, y=MDS2)) + geom_point(aes_string(color="`sample_environment-feature`")) + theme_bw() +
  scale_x_continuous(limits=c(-0.1,0.1)) 


## ---- differential_taxa--------------------------------------------------
library(reshape2)
library(dplyr)
library(MASS)
remerged <- cbind(goslim[,c(1,2,3)], as.data.frame(normalized))

longform <- melt(remerged, id.vars=c("accession", "description", "category"), variable.name="sample", value.name="abund")

longform <- merge(longform, full_metadata[,c("analysis_accession","sample_environment-feature")], by.x="sample", by.y="analysis_accession")

#r <- longform  %>%
#  filter(`sample_environment-feature` %in% c("hydrothermal vent", "sea shore")) %>%
#  mutate(transformed_abund = abund+1) %>%
#  group_by(accession) %>%
#  do(
#    broom::tidy(glm.nb(transformed_abund ~ `sample_environment-feature`, .))
#    ) 




## ---- get_download_urls, results='hide',message=FALSE--------------------
#Find list of available downloads, and filter for 
dl_urls <- mgnify_get_download_urls(mg, full_metadata$analysis_accession, accession_type = "analyses")

## ----show_tgt_urls-------------------------------------------------------
target_urls <- dl_urls[dl_urls$attributes.description.label == "Predicted CDS with annotation",]
head(target_urls)

## ---- list_descriptions--------------------------------------------------
table(dl_urls$attributes.description.label)

## ----get_genome_urls, results='hide'-------------------------------------
genome_urls <- mgnify_get_download_urls(mg, "MGYG-HGUT-04644", accession_type = "genomes")

## ----show_genome_urls----------------------------------------------------
genome_urls[c("id","attributes.file.format.name","download_url")]

## ---- filter_dl_urls, echo=T, message=FALSE------------------------------
#Just select a single file from the target_urls list for demonstration.

#Default behaviour - use local cache.
cached_location = mgnify_download(mg, target_urls$download_url[[41]])

#Specifying a target_filename
specified_location = mgnify_download(mg, target_urls$download_url[[41]], target_filename = "Coding_Sequences_1.fasta.gz")

#Where are the files?
c(cached_location,specified_location)

## ---- simple_parse_function----------------------------------------------
library(Biostrings)

#Simple function to a count of unique sequences matching PFAM amoC/mmoC motif
getAmoCseqs <- function(fname){
  sequences <- readAAStringSet(fname)
  tgtvec <- grepl("PF04896", names(sequences))
  as.data.frame(as.list(table(as.character(sequences[tgtvec]))))
}

## ----do_download_with_read-----------------------------------------------

#Just download a single accession for demonstration, specifying a read_function
amoC_seq_counts <- mgnify_download(mg, target_urls$download_url[[41]] , read_func = getAmoCseqs, usecache = F)

amoC_seq_counts %>% t  


