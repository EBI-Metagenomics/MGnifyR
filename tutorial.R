#restart R
# .rs.restartR()

#clear terminals screen
cat("\014")

#clear variables
rm(list=ls())

# 2.2.1 Assay data 
library(mia)
data(GlobalPatterns, package="mia")
tse <- GlobalPatterns
assays(tse)
assay(tse, "counts")
assay(tse, "counts")[1:5,1:7]

tse <- relAbundanceCounts(tse)
assays(tse)
assay(tse, "relabundance")[1:5,1:7]

# 2.2.2 colData (contains data on the samples)
colData(tse)

# 2.2.3 rowDara (contains data on the features of the analyzed samples)
rowData(tse)

# 2.2.4 rowTree (able to keep track of feature and node relations via two functions)
rowTree(tse)
rowLinks(tse)

# 2.2.5 Alternative experiments
tse_phylum <- agglomerateByRank(tse, "Phylum")
dim(tse)
dim(tse_phylum)
altExp(tse, "Phylum") <- tse_phylum
altExpNames(tse)
tse[,1:10]
dim(altExp(tse[,1:10],"Phylum"))

#2.4.1 Package data
library(mia)
#data(package="mia") #open a window
data("GlobalPatterns", package="mia")
GlobalPatterns

# 2.4.2 ExperimentHub data
library(microbiomeDataSets)
availableDataSets()
mae <- HintikkaXOData


# 4.1.1 Tidy data
library(mia)
data(GlobalPatterns, package="mia")
tse <- GlobalPatterns
tse <- transformSamples(tse, method="relabundance")

molten_tse <- meltAssay(tse,
                        add_row_data = TRUE,
                        add_col_data = TRUE,
                        abund_values = "relabundance")
molten_tse

# 4.1.2 Subsetting
data("GlobalPatterns", package="mia")
tse <- GlobalPatterns
dim(tse)

# 4.1.2.1 Subset by sample (column-wise)
unique(tse$SampleType) # inspect possible values for SampleType
library(magrittr)
library(dplyr)  
tse$SampleType %>% table() # show recurrence for each value
tse_subset_by_sample <- tse[ , tse$SampleType %in% c("Feces", "Skin", "Tongue")]
dim(tse_subset_by_sample)

# 4.1.2.2 Subset by feature (row-wise)
unique(rowData(tse)$Phylum)
rowData(tse)$Phylum %>% tablee()
# 4.1.2.2.1 Non-agglomerated data
tse_subset_by_feature <- tse[rowData(tse)$Phylum %in% c("Actinobacteria", "Chlamydiae") & !is.na(rowData(tse)$Phylum), ]
dim(tse_subset_by_feature)
#4.1.2.2.2 agglomerated data
tse_phylum <- tse %>% agglomerateByRank(rank = "Phylum")
tse_phylum_subset_by_feature <- tse_phylum[rowData(tse_phylum)$Phylum %in% c("Actinobacteria", "Chlamydiae") & !is.na(rowData(tse_phylum)$Phylum), ]
dim(tse_phylum_subset_by_feature)

# 4.1.3 Subset by sample and feature
tse_subset_by_sample_feature <- tse[rowData(tse)$Phylum %in% c("Actinobacteria", "Chlamydiae") & !is.na(rowData(tse)$Phylum), tse$SampleType %in% c("Feces", "Skin", "Tongue")]
dim(tse_subset_by_sample_feature)

# 4.1.2.4 Remove empty columns and rows
tse_genus <- agglomerateByRank(tse, rank = "Genus")
  # List bacteria that we want to include
genera <- c("Class:Thermoprotei", "Genus:Sulfolobus", "Genus:Sediminicola")
tse_genus_sub <- tse_genus[genera, ]
tse_genus_sub
  # List total counts of each sample
colSums(assay(tse_genus_sub, "counts"))
  # Remove samples that do not have present any bacteria
tse_genus_sub <- tse_genus_sub[ , colSums(assay(tse_genus_sub, "counts")) != 0 ]
tse_genus_sub
  # Take only those samples that are collected from feces, skin, or tongue
tse_genus_sub <- tse_genus[ , colData(tse_genus)$SampleType %in% c("Feces", "Skin", "Tongue")]
tse_genus_sub
  # How many bacteria there are that are not present?
sum(rowSums(assay(tse_genus_sub, "counts")) == 0)
  # Take only those bacteria that are present
tse_genus_sub <- tse_genus_sub[rowSums(assay(tse_genus_sub, "counts")) > 0, ]
tse_genus_sub


# 9.1 Visualizing taxonomic composition
library(mia)
data("GlobalPatterns", package="mia")
tse <- GlobalPatterns

# 9.1.1 Composition barplot
library(miaViz)
  # Computing relative abundance
tse <- relAbundanceCounts(tse)
  # Getting top taxa on a Phylum level
tse_phylum <- agglomerateByRank(tse, rank ="Phylum", onRankOnly=TRUE)
top_taxa <- getTopTaxa(tse_phylum,top = 5, abund_values = "relabundance")
  # Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
phylum_renamed <- lapply(rowData(tse)$Phylum, function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(tse)$Phylum <- as.character(phylum_renamed)
  # Visualizing the composition barplot, with samples order by "Bacteroidetes"
plotAbundance(tse, abund_values="relabundance", rank = "Phylum",order_rank_by="abund", order_sample_by = "Bacteroidetes")

# 9.1.2 Composition heatmap
library(ggplot2)
# Add clr-transformation on samples
tse_phylum <- transformSamples(tse_phylum, method = "clr", pseudocount = 1)
# Add z-transformation on features (taxa)
tse_phylum <- transformFeatures(tse_phylum, abund_values = "clr", 
                                method = "z", name = "clr_z")
  # Melts the assay
df <- meltAssay(tse_phylum, abund_values = "clr_z")
  # Determines the scaling of colours
maxval <- round(max(abs(df$clr_z)))
limits <- c(-maxval, maxval)
breaks <- seq(from = min(limits), to = max(limits), by = 0.5)
colours <- c("darkblue", "blue", "white", "red", "darkred")
  # Creates a ggplot object
ggplot(df, aes(x = SampleID, y = FeatureID, fill = clr_z)) +
  geom_tile() +
  scale_fill_gradientn(name = "CLR + Z transform", 
                       breaks = breaks, limits = limits, colours = colours) + 
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.key.size = unit(1, "cm")) +
  labs(x = "Samples", y = "Taxa")


# 13 Multy-assay analyses
library(mia)
  # Load the data
mae <- microbiomeDataSets::HintikkaXOData()
mae
if(!require(stringr)){
  install.packages("stringr")
  library(stringr)
}
  # Drop off those bacteria that do not include information in Phylum or lower levels
mae[[1]] <- mae[[1]][!is.na(rowData(mae[[1]])$Phylum), ]
  # Clean taxonomy data, so that names do not include additional characters
rowData(mae[[1]]) <- DataFrame(apply(rowData(mae[[1]]), 2, 
                                     str_remove, pattern = "._[0-9]__"))
  # Microbiome data
mae[[1]]
  # Metabolite data
mae[[2]]
  # Biomarker data
mae[[3]]
  # Agglomerate microbiome data at Genus level
mae[[1]] <- agglomerateByPrevalence(mae[[1]], rank = "Genus")
  # Does log10 transform for microbiome data
mae[[1]] <- transformSamples(mae[[1]], method = "log10", pseudocount = 1)
  # Give unique names so that we do not have problems when we are creating a plot
rownames(mae[[1]]) <- getTaxonomyLabels(mae[[1]])
  # Cross correlates data sets
correlations <- testExperimentCrossCorrelation(mae, 
                                               experiment1 = 1,
                                               experiment2 = 2,
                                               abund_values1 = "log10", 
                                               abund_values2 = "nmr",
                                               method = "spearman", 
                                               p_adj_threshold = 0.05,
                                               cor_threshold = 0,
                                               # Remove when mia is fixed
                                               # mode = "matrix",
                                               # sort = TRUE,
                                               show_warnings = FALSE)
# if( !require("ComplexHeatmap") ){
#     BiocManager::install("ComplexHeatmap")
#     library("ComplexHeatmap")
# }
#
# # Create a heatmap and store it
# plot <- Heatmap(correlations$cor, 
#                 # Print values to cells
#                 cell_fun = function(j, i, x, y, width, height, fill) {
#                     # If the p-value is under threshold
#                     if( correlations$p_adj[i, j] < 0.05 )
#                         # Print "X"
#                         grid.text(sprintf("%s", "X"), x, y, gp = gpar(fontsize = 8, col = "black"))
#                     },
#                 heatmap_legend_param = list(title = "", legend_height = unit(5, "cm"))
#                 )
# plot
library(ggplot2)

  # Determines the scaling of colors
limits <- c(-1, 1)
breaks <- seq(from = min(limits), to = max(limits), by = 0.2)
colours <- c("darkblue", "blue", "white", "red", "darkred")

  # Which observation have p-value under 0.05? --> creates a subset
cor_table_sub <- correlations[which(correlations[["p_adj"]] < 0.05), ]

  # Creates a ggplot object
ggplot(correlations, aes(x = Var1, y = Var2, fill = cor)) +
  geom_tile() +
  
  scale_fill_gradientn(name = "Correlation",
                       breaks = breaks, limits = limits, colours = colours) +
  
  # Adds label to those observations that have p-value under 0.05
  geom_text(data = cor_table_sub, aes(x = Var1, y = Var2, label = "+")) +
  
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.key.size = unit(1, "cm")) +
  labs(x = "Taxa", y = "Metabolites")
