library(vegan)
library(ggplot2)
library(phyloseq)

library(MGnifyR)

mg <- mgnify_client(usecache = T, cache_dir = '/tmp/MGnify_demo')


####### Queries:
mgnify_query(mg, "studies", biome_name="Wastewater", usecache = T)
mgnify_query(mg, "samples", latitude_gte=60.0, experiment_type="metagenomic", usecache = T)
m <- mgnify_query(mg, "samples", study_accession = "MGYS00003725", usecache=T)
acc_list <- mgnify_analyses_from_samples(mg, m$accession)
df <- mgnify_get_analyses_metadata(mg, acc_list)
df


##### Single study retrieval
#Amplicon: Oral health of young adults: Amplicon study
om_analyses <- mgnify_analyses_from_studies(mg, "MGYS00002277")
om_metadata_df <- mgnify_get_analyses_metadata(mg, om_analyses)
t(head(om_metadata_df))

om_ps <- mgnify_get_analyses_phyloseq(mg, om_analyses, tax_SU = "SSU")

om_ps_sub <- subset_samples(om_ps, sample_sums(om_ps) > 10000)

omps <- rarefy_even_depth(om_ps_sub)
omps

#plt1 <- plot_bar(omps, fill="Class", facet_grid = "sample_sample.desc") + theme(legend.position = "none")

alpha_div <- estimate_richness(omps)

adf <- cbind.data.frame(sample_data(omps)$`sample_sample.desc`, alpha_div$InvSimpson)
colnames(adf) <- c("factor","value")

ggplot(adf, aes(x=factor, y=value)) + geom_boxplot(width=0.1) + geom_jitter(width=0.1) + theme_bw()

omps_ord <- ordinate(omps, method = "PCoA" , distance = "bray")
plot_ordination(omps, omps_ord, color = "sample_sample.desc") + theme_bw()



##### Multi-biome metagenome analysis
set.seed(11)
mg <- mgnify_client(usecache = T, cache_dir = "/tmp/mgnify_cache")

#Study accessions

#Saltmarsh metagenomes HiSeq4000 : MGYS00001447 - 48 samples
#Healthy human gut metagenomes : MGYS00001442 30 odd samples
#Marine Subseafloor microbes at Mid-Cayman Rise: MGYS00001282


soil <- mgnify_analyses_from_studies(mg, "MGYS00001447")
human <- mgnify_analyses_from_studies(mg, "MGYS00001442")
seafloor <- mgnify_analyses_from_studies(mg, "MGYS00001282")
seafloor <- sample(seafloor, 40)

accessions <- c(soil,human,seafloor)

metadata <- mgnify_get_analyses_metadata(mg, accessions)
head(metadata)

goterms <-  mgnify_get_analyses_results(mg, accessions, retrievelist = "go-slim")$`go-slim`


m <- goterms[,c(-1,-2,-3)]

normed_m <- apply(m, 2, function(x) x/sum(x))

nmds <- vegan::metaMDS(t(normed_m))

pltdat <- as.data.frame(scores(vare.mds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
pltdat$grp <- metadata[rownames(pltdat),"study_accession"]  #  add the grp variable created earlier

ggplot() + geom_point(data=pltdat,aes(x=NMDS1,y=NMDS2,colour=grp),size=3) + theme_bw()




