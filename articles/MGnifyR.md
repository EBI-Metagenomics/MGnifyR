# MGnifyR

## Introduction

`MGnifyR` fetches microbiome data from
[MGnify](https://www.ebi.ac.uk/metagenomics), EBI’s resource for the
analysis of metagenomic, metatranscriptomic and amplicon data. In
MGnify,

- a **study** (e.g., `MGYS00005116`) is a collection of samples,
- a **sample** (e.g., `SAMN02595037`) is a biological specimen, e.g., a
  stool or soil sample,
- a **run** (e.g., `SRR1144937`) is a sequencing of a sample, i.e., its
  reads. A sample can have multiple runs, e.g., an amplicon and a
  metagenome.
- an **assembly** (e.g., `ERZ1746131`) is made of the reads of one or
  more runs, and
- an **analysis** (e.g., `MGYA00383350`) is a run or an assembly
  processed with an MGnify **pipeline**. The data is fetched by
  analysis, so each column of the result is an analysis. A run can be
  analysed with multiple pipeline versions, so it can have multiple
  analyses. An analysis usually has one run, but some have multiple
  runs, e.g., runs that were analysed together.

`MGnifyR` returns the data as a `TreeSummarizedExperiment`, so that it
can be analysed directly with the [Orchestrating Microbiome Analysis
with Bioconductor (OMA)](https://microbiome.github.io/) tools. See the
[OMA book](https://microbiome.github.io/OMA/) for analysis examples.

## Installation

\
`BiocManager``::`[`install`](https://bioconductor.github.io/BiocManager/reference/install.html)`(``"MGnifyR"``)`

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`MGnifyR`](https://github.com/EBI-Metagenomics/MGnifyR)`)`\
`#> Loading required package: TreeSummarizedExperiment`\
`#> Loading required package: SingleCellExperiment`\
`#> Loading required package: SummarizedExperiment`\
`#> Loading required package: MatrixGenerics`\
`#> Loading required package: matrixStats`\
`#> `\
`#> Attaching package: 'MatrixGenerics'`\
`#> The following objects are masked from 'package:matrixStats':`\
`#> `\
`#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,`\
`#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,`\
`#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,`\
`#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,`\
`#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,`\
`#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,`\
`#>     colWeightedMeans, colWeightedMedians, colWeightedSds,`\
`#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,`\
`#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,`\
`#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,`\
`#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,`\
`#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,`\
`#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,`\
`#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,`\
`#>     rowWeightedSds, rowWeightedVars`\
`#> Loading required package: GenomicRanges`\
`#> Loading required package: stats4`\
`#> Loading required package: BiocGenerics`\
`#> Loading required package: generics`\
`#> `\
`#> Attaching package: 'generics'`\
`#> The following objects are masked from 'package:base':`\
`#> `\
`#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,`\
`#>     setequal, union`\
`#> `\
`#> Attaching package: 'BiocGenerics'`\
`#> The following objects are masked from 'package:stats':`\
`#> `\
`#>     IQR, mad, sd, var, xtabs`\
`#> The following object is masked from 'package:utils':`\
`#> `\
`#>     data`\
`#> The following objects are masked from 'package:base':`\
`#> `\
`#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,`\
`#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,`\
`#>     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,`\
`#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,`\
`#>     rbind, Reduce, rownames, sapply, saveRDS, scale, sequence, table,`\
`#>     tapply, transform, unique, unsplit, which.max, which.min`\
`#> Loading required package: S4Vectors`\
`#> `\
`#> Attaching package: 'S4Vectors'`\
`#> The following object is masked from 'package:utils':`\
`#> `\
`#>     findMatches`\
`#> The following objects are masked from 'package:base':`\
`#> `\
`#>     expand.grid, I, unname`\
`#> Loading required package: IRanges`\
`#> Loading required package: Seqinfo`\
`#> Loading required package: Biobase`\
`#> Welcome to Bioconductor`\
`#> `\
`#>     Vignettes contain introductory material; view with`\
`#>     'browseVignettes()'. To cite Bioconductor, see`\
`#>     'citation("Biobase")', and for packages 'citation("pkgname")'.`\
`#> `\
`#> Attaching package: 'Biobase'`\
`#> The following object is masked from 'package:MatrixGenerics':`\
`#> `\
`#>     rowMedians`\
`#> The following objects are masked from 'package:matrixStats':`\
`#> `\
`#>     anyMissing, rowMedians`\
`#> Loading required package: Biostrings`\
`#> Loading required package: XVector`\
`#> `\
`#> Attaching package: 'Biostrings'`\
`#> The following object is masked from 'package:base':`\
`#> `\
`#>     strsplit`

## Browse studies

[`searchStudies()`](../reference/searchStudies.md) returns a table of
studies. Studies can be filtered by biome, pipeline version and text
search, and sorted with `order`. The available biomes are listed in the
[MGnify website](https://www.ebi.ac.uk/metagenomics/browse/biomes).

\
`studies`` ``<-`` `[`searchStudies`](../reference/searchStudies.md)`(`\
`    biome ``=`` ``"root:Engineered:Built environment"``, max.hits ``=`` ``10``)`\
`studies``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"accession"``, ``"title"``, ``"biome_lineage"``)``]`\
`#>       accession`\
`#> 1  MGYS00001928`\
`#> 2  MGYS00002607`\
`#> 3  MGYS00004734`\
`#> 4  MGYS00005116`\
`#> 5  MGYS00005117`\
`#> 6  MGYS00005289`\
`#> 7  MGYS00005290`\
`#> 8  MGYS00005291`\
`#> 9  MGYS00005292`\
`#> 10 MGYS00005293`\
`#>                                                                                                                                   title`\
`#> 1                                                                                                         Hot water plumbing microbiome`\
`#> 2                                                                                                          Marine Aggregates Metagenome`\
`#> 3                                                                                                                       ISS Metagenomes`\
`#> 4                                                                  International Space Station vacuum cleaner bag Targeted Locus (Loci)`\
`#> 5                                                                             Analysis of dust samples from the Russian part of the ISS`\
`#> 6  Preparing for the manned Mars journey: Microbiome dynamics in the confined MARS 500 habitat during simulated Mars flight and landing`\
`#> 7                                                                               Microbial diversity on kitchen sinks and bathroom tiles`\
`#> 8                                                                                         Inflated Lunar/Mars Analog Habitat Microbiome`\
`#> 9                                                               International Space Station Microbial Observatory - Microbial Diversity`\
`#> 10                                                                            International Space Station flight project EXTREMOPHILES.`\
`#>                        biome_lineage`\
`#> 1  root:Engineered:Built environment`\
`#> 2  root:Engineered:Built environment`\
`#> 3  root:Engineered:Built environment`\
`#> 4  root:Engineered:Built environment`\
`#> 5  root:Engineered:Built environment`\
`#> 6  root:Engineered:Built environment`\
`#> 7  root:Engineered:Built environment`\
`#> 8  root:Engineered:Built environment`\
`#> 9  root:Engineered:Built environment`\
`#> 10 root:Engineered:Built environment`

## Fetch data

Give the studies, or their accessions, to
[`fetchStudies()`](../reference/fetchStudies.md). The taxonomic profiles
of all the studies are fetched and combined into a single
`TreeSummarizedExperiment`. Taxa are matched across studies by their
taxonomy.

\
`tse`` ``<-`` `[`fetchStudies`](../reference/fetchStudies.md)`(``"MGYS00005116"``)`\
`#> Fetching 2 analyses one by one, because they are missing from the summary tables of studies. This can take a while.`\
`#> The data includes analyses from multiple pipeline versions (V4.1, V5, V6). Their profiles are based on different reference databases. Consider subsetting the result of searchAnalyses() and fetching the analyses with fetchAnalyses().`\
`tse`\
`#> class: TreeSummarizedExperiment `\
`#> dim: 1014 6 `\
`#> metadata(0):`\
`#> assays(1): counts`\
`#> rownames(1014): sk__Bacteria sk__Bacteria;k__;p__Acidobacteria ...`\
`#>   sk__Eukaryota;k__Viridiplantae;p__Streptophyta;c__Magnoliopsida;o__Poales;f__Poaceae;g__Triticum`\
`#>   sk__Eukaryota;k__Viridiplantae;p__Streptophyta;c__Magnoliopsida;o__Solanales;f__Solanaceae;g__Nicotiana;s__Nicotiana_attenuata`\
`#> rowData names(8): superkingdom kingdom ... genus species`\
`#> colnames(6): MGYA00383350 MGYA00383351 ... MGYA01026512 MGYA01026513`\
`#> colData names(16): analysis_accession study_accession ...`\
`#>   sample_updated_at abundance_table`\
`#> reducedDimNames(0):`\
`#> mainExpName: NULL`\
`#> altExpNames(0):`\
`#> rowLinks: NULL`\
`#> rowTree: NULL`\
`#> colLinks: NULL`\
`#> colTree: NULL`

The `colData` includes the accessions of the analysis, study, sample and
run, the pipeline version and the other fields of the analyses in
MGnify, e.g., the title of the sample (`sample_sample_title`). More
metadata can be added from other sources: the samples in MGnify
(`samples = TRUE`) and in BioSamples (`biosamples = TRUE`), and the
studies (`studies = TRUE`). BioSamples has all the attributes of the
samples, including study-specific variables, e.g., disease status or
treatment. Each source needs its own requests, so they are not added by
default. The data is cached, so adding a source later fetches only the
new data.

\
`tse`` ``<-`` `[`fetchStudies`](../reference/fetchStudies.md)`(``"MGYS00005116"``, biosamples ``=`` ``TRUE``, studies ``=`` ``TRUE``)`\
`#> Fetching 2 analyses one by one, because they are missing from the summary tables of studies. This can take a while.`\
`#> The data includes analyses from multiple pipeline versions (V4.1, V5, V6). Their profiles are based on different reference databases. Consider subsetting the result of searchAnalyses() and fetching the analyses with fetchAnalyses().`\
`colData``(``tse``)``[``, `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    ``"study_accession"``, ``"pipeline_version"``, ``"sample_accession"``, ``"run_accession"``,`\
`    ``"study_biome_lineage"``)``]`\
`#> DataFrame with 6 rows and 5 columns`\
`#>              study_accession pipeline_version sample_accession run_accession`\
`#>                  <character>      <character>      <character>   <character>`\
`#> MGYA00383350    MGYS00005116             V4.1     SAMN02594912    SRR1139729`\
`#> MGYA00383351    MGYS00005116             V4.1     SAMN02595037    SRR1144937`\
`#> MGYA00571587    MGYS00005116               V5     SAMN02594912    SRR1139729`\
`#> MGYA00571730    MGYS00005116               V5     SAMN02595037    SRR1144937`\
`#> MGYA01026512    MGYS00005116               V6     SAMN02595037    SRR1144937`\
`#> MGYA01026513    MGYS00005116               V6     SAMN02594912    SRR1139729`\
`#>                 study_biome_lineage`\
`#>                         <character>`\
`#> MGYA00383350 root:Engineered:Buil..`\
`#> MGYA00383351 root:Engineered:Buil..`\
`#> MGYA00571587 root:Engineered:Buil..`\
`#> MGYA00571730 root:Engineered:Buil..`\
`#> MGYA01026512 root:Engineered:Buil..`\
`#> MGYA01026513 root:Engineered:Buil..`\
`# Columns of each source`\
[`table`](https://rdrr.io/r/base/table.html)`(`[`sub`](https://rdrr.io/r/base/grep.html)`(``"_.*"``, ``""``, `[`colnames`](https://rdrr.io/r/base/colnames.html)`(``colData``(``tse``)``)``)``)`\
`#> `\
`#>  abundance   analysis   assembly biosamples experiment instrument   pipeline `\
`#>          1          1          1         24          1          2          1 `\
`#>        run     sample      study `\
`#>          4          4         13`

The `rowData` includes the taxonomy.

\
[`head`](https://rdrr.io/r/utils/head.html)`(``rowData``(``tse``)``)`\
`#> DataFrame with 6 rows and 8 columns`\
`#>                                                                                                                superkingdom`\
`#>                                                                                                                 <character>`\
`#> sk__Bacteria                                                                                                       Bacteria`\
`#> sk__Bacteria;k__;p__Acidobacteria                                                                                  Bacteria`\
`#> sk__Bacteria;k__;p__Acidobacteria;c__Acidobacteriia;o__Acidobacteriales;f__Acidobacteriaceae                       Bacteria`\
`#> sk__Bacteria;k__;p__Acidobacteria;c__Holophagae;o__Holophagales                                                    Bacteria`\
`#> sk__Bacteria;k__;p__Actinobacteria                                                                                 Bacteria`\
`#> sk__Bacteria;k__;p__Actinobacteria;c__Acidimicrobiia;o__Acidimicrobiales;f__Acidimicrobiaceae;g__Ilumatobacter     Bacteria`\
`#>                                                                                                                    kingdom`\
`#>                                                                                                                <character>`\
`#> sk__Bacteria                                                                                                            NA`\
`#> sk__Bacteria;k__;p__Acidobacteria                                                                                       NA`\
`#> sk__Bacteria;k__;p__Acidobacteria;c__Acidobacteriia;o__Acidobacteriales;f__Acidobacteriaceae                            NA`\
`#> sk__Bacteria;k__;p__Acidobacteria;c__Holophagae;o__Holophagales                                                         NA`\
`#> sk__Bacteria;k__;p__Actinobacteria                                                                                      NA`\
`#> sk__Bacteria;k__;p__Actinobacteria;c__Acidimicrobiia;o__Acidimicrobiales;f__Acidimicrobiaceae;g__Ilumatobacter          NA`\
`#>                                                                                                                        phylum`\
`#>                                                                                                                   <character>`\
`#> sk__Bacteria                                                                                                               NA`\
`#> sk__Bacteria;k__;p__Acidobacteria                                                                               Acidobacteria`\
`#> sk__Bacteria;k__;p__Acidobacteria;c__Acidobacteriia;o__Acidobacteriales;f__Acidobacteriaceae                    Acidobacteria`\
`#> sk__Bacteria;k__;p__Acidobacteria;c__Holophagae;o__Holophagales                                                 Acidobacteria`\
`#> sk__Bacteria;k__;p__Actinobacteria                                                                             Actinobacteria`\
`#> sk__Bacteria;k__;p__Actinobacteria;c__Acidimicrobiia;o__Acidimicrobiales;f__Acidimicrobiaceae;g__Ilumatobacter Actinobacteria`\
`#>                                                                                                                         class`\
`#>                                                                                                                   <character>`\
`#> sk__Bacteria                                                                                                               NA`\
`#> sk__Bacteria;k__;p__Acidobacteria                                                                                          NA`\
`#> sk__Bacteria;k__;p__Acidobacteria;c__Acidobacteriia;o__Acidobacteriales;f__Acidobacteriaceae                   Acidobacteriia`\
`#> sk__Bacteria;k__;p__Acidobacteria;c__Holophagae;o__Holophagales                                                    Holophagae`\
`#> sk__Bacteria;k__;p__Actinobacteria                                                                                         NA`\
`#> sk__Bacteria;k__;p__Actinobacteria;c__Acidimicrobiia;o__Acidimicrobiales;f__Acidimicrobiaceae;g__Ilumatobacter Acidimicrobiia`\
`#>                                                                                                                           order`\
`#>                                                                                                                     <character>`\
`#> sk__Bacteria                                                                                                                 NA`\
`#> sk__Bacteria;k__;p__Acidobacteria                                                                                            NA`\
`#> sk__Bacteria;k__;p__Acidobacteria;c__Acidobacteriia;o__Acidobacteriales;f__Acidobacteriaceae                   Acidobacteriales`\
`#> sk__Bacteria;k__;p__Acidobacteria;c__Holophagae;o__Holophagales                                                    Holophagales`\
`#> sk__Bacteria;k__;p__Actinobacteria                                                                                           NA`\
`#> sk__Bacteria;k__;p__Actinobacteria;c__Acidimicrobiia;o__Acidimicrobiales;f__Acidimicrobiaceae;g__Ilumatobacter Acidimicrobiales`\
`#>                                                                                                                           family`\
`#>                                                                                                                      <character>`\
`#> sk__Bacteria                                                                                                                  NA`\
`#> sk__Bacteria;k__;p__Acidobacteria                                                                                             NA`\
`#> sk__Bacteria;k__;p__Acidobacteria;c__Acidobacteriia;o__Acidobacteriales;f__Acidobacteriaceae                   Acidobacteriaceae`\
`#> sk__Bacteria;k__;p__Acidobacteria;c__Holophagae;o__Holophagales                                                               NA`\
`#> sk__Bacteria;k__;p__Actinobacteria                                                                                            NA`\
`#> sk__Bacteria;k__;p__Actinobacteria;c__Acidimicrobiia;o__Acidimicrobiales;f__Acidimicrobiaceae;g__Ilumatobacter Acidimicrobiaceae`\
`#>                                                                                                                        genus`\
`#>                                                                                                                  <character>`\
`#> sk__Bacteria                                                                                                              NA`\
`#> sk__Bacteria;k__;p__Acidobacteria                                                                                         NA`\
`#> sk__Bacteria;k__;p__Acidobacteria;c__Acidobacteriia;o__Acidobacteriales;f__Acidobacteriaceae                              NA`\
`#> sk__Bacteria;k__;p__Acidobacteria;c__Holophagae;o__Holophagales                                                           NA`\
`#> sk__Bacteria;k__;p__Actinobacteria                                                                                        NA`\
`#> sk__Bacteria;k__;p__Actinobacteria;c__Acidimicrobiia;o__Acidimicrobiales;f__Acidimicrobiaceae;g__Ilumatobacter Ilumatobacter`\
`#>                                                                                                                    species`\
`#>                                                                                                                <character>`\
`#> sk__Bacteria                                                                                                            NA`\
`#> sk__Bacteria;k__;p__Acidobacteria                                                                                       NA`\
`#> sk__Bacteria;k__;p__Acidobacteria;c__Acidobacteriia;o__Acidobacteriales;f__Acidobacteriaceae                            NA`\
`#> sk__Bacteria;k__;p__Acidobacteria;c__Holophagae;o__Holophagales                                                         NA`\
`#> sk__Bacteria;k__;p__Actinobacteria                                                                                      NA`\
`#> sk__Bacteria;k__;p__Actinobacteria;c__Acidimicrobiia;o__Acidimicrobiales;f__Acidimicrobiaceae;g__Ilumatobacter          NA`

### Pipeline versions

The taxonomies of different pipeline versions are based on different
reference databases, so they might not be directly comparable. To fetch
analyses of a single pipeline version, select them with
[`searchAnalyses()`](../reference/searchStudies.md) as shown below.

## Types of profiles

By default, small subunit rRNA (`"SSU"`) profiles are fetched. `type`
tells which profile is fetched, e.g., `"mOTUs"` for metagenomes and
metatranscriptomes, large subunit rRNA (`"LSU"`), `"UNITE"`,
`"ITSoneDB"` or `"PR2"`. Only the analyses that have the type are
fetched. See [`?fetchStudies`](../reference/fetchStudies.md) for all the
types.

\
`tse_lsu`` ``<-`` `[`fetchStudies`](../reference/fetchStudies.md)`(``"MGYS00005116"``, type ``=`` ``"LSU"``)`\
`#> Fetching 5 analyses one by one, because they are missing from the summary tables of studies. This can take a while.`\
`#> No LSU abundance table could be fetched for 3 analyses (no such table, or it could not be downloaded): 'MGYA00383350', 'MGYA00383351', 'MGYA00571587'`\
`#> The data includes analyses from multiple pipeline versions (V5, V6). Their profiles are based on different reference databases. Consider subsetting the result of searchAnalyses() and fetching the analyses with fetchAnalyses().`\
`tse_lsu`\
`#> class: TreeSummarizedExperiment `\
`#> dim: 2 3 `\
`#> metadata(0):`\
`#> assays(1): counts`\
`#> rownames(2):`\
`#>   sk__Eukaryota;k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales`\
`#>   Unclassified`\
`#> rowData names(8): superkingdom kingdom ... genus species`\
`#> colnames(3): MGYA00571730 MGYA01026512 MGYA01026513`\
`#> colData names(16): analysis_accession study_accession ...`\
`#>   sample_updated_at abundance_table`\
`#> reducedDimNames(0):`\
`#> mainExpName: NULL`\
`#> altExpNames(0):`\
`#> rowLinks: NULL`\
`#> rowTree: NULL`\
`#> colLinks: NULL`\
`#> colTree: NULL`

## Functional profiles

Metagenomes, metatranscriptomes and assemblies also have functional
profiles:

- `"InterPro"`, `"GO"` and `"GO-slim"`: pipelines v1-v6
- `"Pfam"` and `"KO"` (KEGG orthologs): pipelines v5-v6
- `"antiSMASH"` (biosynthetic gene clusters): assemblies of pipelines
  v5-v6

In pipeline v6, only assemblies have functional profiles. Amplicons have
no functional profiles. The rows are, e.g., GO terms, and `rowData` has
their descriptions.

\
`analyses`` ``<-`` `[`searchAnalyses`](../reference/searchStudies.md)`(``"MGYS00001442"``)`\
`analyses`` ``<-`` ``analyses``[``1``:``3``, ``]`\
`tse_go`` ``<-`` `[`fetchAnalyses`](../reference/fetchStudies.md)`(``analyses``, type ``=`` ``"GO-slim"``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(``rowData``(``tse_go``)``)`\
`#> DataFrame with 6 rows and 2 columns`\
`#>                       description           category`\
`#>                       <character>        <character>`\
`#> GO:0000015 phosphopyruvate hydr.. cellular_component`\
`#> GO:0000150   recombinase activity molecular_function`\
`#> GO:0000160 phosphorelay signal .. biological_process`\
`#> GO:0000166     nucleotide binding molecular_function`\
`#> GO:0003674     molecular function molecular_function`\
`#> GO:0003676   nucleic acid binding molecular_function`

## Multiple profiles

To fetch multiple types, give them all to `type`. The result is a
`MultiAssayExperiment` with a `TreeSummarizedExperiment` for each type.
The profiles are linked by the sample that was analysed, so, e.g., the
taxonomic profile of an amplicon and the functional profile of a
metagenome of the same sample are linked. Types that have no data are
left out.

\
`mae`` ``<-`` `[`fetchAnalyses`](../reference/fetchStudies.md)`(``analyses``, type ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"SSU"``, ``"GO-slim"``, ``"InterPro"``)``)`\
`mae`\
`#> A MultiAssayExperiment object of 3 listed`\
`#>  experiments with user-defined names and respective classes.`\
`#>  Containing an ExperimentList class object of length 3:`\
`#>  [1] SSU: TreeSummarizedExperiment with 441 rows and 3 columns`\
`#>  [2] GO-slim: TreeSummarizedExperiment with 116 rows and 3 columns`\
`#>  [3] InterPro: TreeSummarizedExperiment with 8495 rows and 3 columns`\
`#> Functionality:`\
`#>  experiments() - obtain the ExperimentList instance`\
`#>  colData() - the primary/phenotype DataFrame`\
`#>  sampleMap() - the sample coordination DataFrame`\
`` #>  `$`, `[`, `[[` - extract colData columns, subset, or experiment ``\
`#>  *Format() - convert into a long or wide DataFrame`\
`#>  assays() - convert ExperimentList to a SimpleList of matrices`\
`#>  exportClass() - save data to flat files`

Each profile is a `TreeSummarizedExperiment`, and the sample metadata is
in `colData`.

\
`mae``[[``"GO-slim"``]``]`\
`#> class: TreeSummarizedExperiment `\
`#> dim: 116 3 `\
`#> metadata(0):`\
`#> assays(1): counts`\
`#> rownames(116): GO:0000015 GO:0000150 ... GO:1902494 GO:1990204`\
`#> rowData names(2): description category`\
`#> colnames(3): MGYA00096976 MGYA00096977 MGYA00096978`\
`#> colData names(16): analysis_accession study_accession ...`\
`#>   sample_updated_at abundance_table`\
`#> reducedDimNames(0):`\
`#> mainExpName: NULL`\
`#> altExpNames(0):`\
`#> rowLinks: NULL`\
`#> rowTree: NULL`\
`#> colLinks: NULL`\
`#> colTree: NULL`\
`colData``(``mae``)`\
`#> DataFrame with 3 rows and 8 columns`\
`#>               study_accession sample_accession run_experiment_type`\
`#>                   <character>      <character>         <character>`\
`#> SAMEA53400418    MGYS00001442    SAMEA53400418         Metagenomic`\
`#> SAMEA53410918    MGYS00001442    SAMEA53410918         Metagenomic`\
`#> SAMEA53405668    MGYS00001442    SAMEA53405668         Metagenomic`\
`#>               run_sample_accession run_study_accession  sample_ena_accessions`\
`#>                        <character>         <character>            <character>`\
`#> SAMEA53400418        SAMEA53400418        MGYS00001442 ["SAMEA53400418","ER..`\
`#> SAMEA53410918        SAMEA53410918        MGYS00001442 ["SAMEA53410918","ER..`\
`#> SAMEA53405668        SAMEA53405668        MGYS00001442 ["ERS1509322","SAMEA..`\
`#>                sample_sample_title      sample_updated_at`\
`#>                        <character>            <character>`\
`#> SAMEA53400418 human gut metagenome 2026-10-31T12:19:07...`\
`#> SAMEA53410918 human gut metagenome 2026-10-31T12:18:47...`\
`#> SAMEA53405668 human gut metagenome 2026-10-31T12:18:46...`

## Select analyses

To select specific analyses, e.g., based on sample metadata, use
[`searchAnalyses()`](../reference/searchStudies.md). It returns one row
per analysis with the same columns as the `colData`. Filter the table
and give it to [`fetchAnalyses()`](../reference/fetchStudies.md).

\
`analyses`` ``<-`` `[`searchAnalyses`](../reference/searchStudies.md)`(``"MGYS00005116"``)`\
`analyses``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"analysis_accession"``, ``"pipeline_version"``, ``"sample_accession"``)``]`\
`#>   analysis_accession pipeline_version sample_accession`\
`#> 1       MGYA01026512               V6     SAMN02595037`\
`#> 2       MGYA01026513               V6     SAMN02594912`\
`#> 3       MGYA00383350             V4.1     SAMN02594912`\
`#> 4       MGYA00383351             V4.1     SAMN02595037`\
`#> 5       MGYA00571587               V5     SAMN02594912`\
`#> 6       MGYA00571730               V5     SAMN02595037`\
\
`analyses`` ``<-`` ``analyses``[``analyses``[[``"pipeline_version"``]``]`` ``==`` ``"V6"``, ``]`\
`tse`` ``<-`` `[`fetchAnalyses`](../reference/fetchStudies.md)`(``analyses``)`\
`#> Fetching 2 analyses one by one, because they are missing from the summary tables of studies. This can take a while.`\
`tse`\
`#> class: TreeSummarizedExperiment `\
`#> dim: 520 2 `\
`#> metadata(0):`\
`#> assays(1): counts`\
`#> rownames(520):`\
`#>   sk__Bacteria;k__;p__Actinomycetota;c__Actinomycetes;o__Propionibacteriales;f__Propionibacteriaceae;g__Cutibacterium`\
`#>   sk__Bacteria;k__;p__Bacillota;c__Bacilli;o__Bacillales;f__Staphylococcaceae;g__Staphylococcus`\
`#>   ...`\
`#>   sk__Eukaryota;k__Viridiplantae;p__Streptophyta;c__Magnoliopsida;o__Poales;f__Poaceae;g__Triticum`\
`#>   sk__Eukaryota;k__Viridiplantae;p__Streptophyta;c__Magnoliopsida;o__Solanales;f__Solanaceae;g__Nicotiana;s__Nicotiana_attenuata`\
`#> rowData names(8): superkingdom kingdom ... genus species`\
`#> colnames(2): MGYA01026512 MGYA01026513`\
`#> colData names(16): analysis_accession study_accession ...`\
`#>   sample_updated_at abundance_table`\
`#> reducedDimNames(0):`\
`#> mainExpName: NULL`\
`#> altExpNames(0):`\
`#> rowLinks: NULL`\
`#> rowTree: NULL`\
`#> colLinks: NULL`\
`#> colTree: NULL`

A study can have several experiment types, e.g., metatranscriptomes and
assemblies. Select them by `experiment_type`.

\
`analyses`` ``<-`` `[`searchAnalyses`](../reference/searchStudies.md)`(``"MGYS00006572"``)`\
[`table`](https://rdrr.io/r/base/table.html)`(``analyses``[[``"experiment_type"``]``]``)`\
`#> `\
`#>           Assembly Metatranscriptomic `\
`#>                  2                  2`\
\
`analyses`` ``<-`` ``analyses``[``analyses``[[``"experiment_type"``]``]`` ``==`` ``"Metatranscriptomic"``, ``]`\
`mae`` ``<-`` `[`fetchAnalyses`](../reference/fetchStudies.md)`(``analyses``, type ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"SSU"``, ``"GO-slim"``)``)`\
`mae`\
`#> A MultiAssayExperiment object of 2 listed`\
`#>  experiments with user-defined names and respective classes.`\
`#>  Containing an ExperimentList class object of length 2:`\
`#>  [1] SSU: TreeSummarizedExperiment with 403 rows and 2 columns`\
`#>  [2] GO-slim: TreeSummarizedExperiment with 116 rows and 2 columns`\
`#> Functionality:`\
`#>  experiments() - obtain the ExperimentList instance`\
`#>  colData() - the primary/phenotype DataFrame`\
`#>  sampleMap() - the sample coordination DataFrame`\
`` #>  `$`, `[`, `[[` - extract colData columns, subset, or experiment ``\
`#>  *Format() - convert into a long or wide DataFrame`\
`#>  assays() - convert ExperimentList to a SimpleList of matrices`\
`#>  exportClass() - save data to flat files`

[`fetchAnalyses()`](../reference/fetchStudies.md) also accepts analysis
accessions, and [`searchSamples()`](../reference/searchStudies.md)
returns the samples of studies.

The arguments `samples`, `biosamples` and `studies` (see Fetch data)
choose the sources of metadata in
[`searchAnalyses()`](../reference/searchStudies.md) too. Also the
details of the analyses (`analyses = TRUE`), e.g., the quality control
of pipeline v6, and of the assemblies (`assemblies = TRUE`) can be
added. They need one request per analysis or assembly.

\
`analyses`` ``<-`` `[`searchAnalyses`](../reference/searchStudies.md)`(``"MGYS00005116"``, analyses ``=`` ``TRUE``)`\
[`table`](https://rdrr.io/r/base/table.html)`(`[`sub`](https://rdrr.io/r/base/grep.html)`(``"_.*"``, ``""``, `[`colnames`](https://rdrr.io/r/base/colnames.html)`(``analyses``)``)``)`\
`#> `\
`#>   analysis   assembly experiment instrument   pipeline        run     sample `\
`#>         48          1          1          2          1          4          4 `\
`#>      study `\
`#>          1`

## Caching and private data

Downloaded data is cached for the R session, so repeated calls are fast.
To keep the cache between sessions, set a cache directory.

\
[`options`](https://rdrr.io/r/base/options.html)`(``MGnifyR.cache ``=`` ``"path/to/cache"``)`

Most MGnify data is public. To access private data, log in with ENA
Webin credentials.

\
[`mgnifyLogin`](../reference/mgnifyLogin.md)`(``"Webin-1122334"``, ``"SecretPassword"``)`

## Session info

\
[`sessionInfo`](https://rdrr.io/r/utils/sessionInfo.html)`(``)`\
`#> R version 4.6.1 (2026-06-24)`\
`#> Platform: x86_64-pc-linux-gnu`\
`#> Running under: Ubuntu 24.04.5 LTS`\
`#> `\
`#> Matrix products: default`\
`#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 `\
`#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0`\
`#> `\
`#> locale:`\
`#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              `\
`#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    `\
`#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   `\
`#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 `\
`#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            `\
`#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       `\
`#> `\
`#> time zone: UTC`\
`#> tzcode source: system (glibc)`\
`#> `\
`#> attached base packages:`\
`#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  `\
`#> [8] base     `\
`#> `\
`#> other attached packages:`\
`#>  [1] MGnifyR_1.99.0                  TreeSummarizedExperiment_2.21.0`\
`#>  [3] Biostrings_2.81.9               XVector_0.53.0                 `\
`#>  [5] SingleCellExperiment_1.35.2     SummarizedExperiment_1.43.0    `\
`#>  [7] Biobase_2.73.2                  GenomicRanges_1.65.4           `\
`#>  [9] Seqinfo_1.3.2                   IRanges_2.47.5                 `\
`#> [11] S4Vectors_0.51.10               BiocGenerics_0.59.12           `\
`#> [13] generics_0.1.4                  MatrixGenerics_1.25.0          `\
`#> [15] matrixStats_1.5.0              `\
`#> `\
`#> loaded via a namespace (and not attached):`\
`#>  [1] yulab.utils_0.2.5           rappdirs_0.3.4             `\
`#>  [3] tidyr_1.3.2                 SparseArray_1.13.4         `\
`#>  [5] lattice_0.23-1              digest_0.6.39              `\
`#>  [7] magrittr_2.0.5              evaluate_1.0.5             `\
`#>  [9] grid_4.6.1                  fastmap_1.2.0              `\
`#> [11] jsonlite_2.0.0              Matrix_1.7-6               `\
`#> [13] ape_5.8-1                   purrr_1.2.2                `\
`#> [15] lazyeval_0.2.3              codetools_0.2-20           `\
`#> [17] abind_1.4-8                 cli_3.6.6                  `\
`#> [19] rlang_1.3.0                 crayon_1.5.3               `\
`#> [21] tidytree_0.4.8              DelayedArray_0.39.8        `\
`#> [23] yaml_2.3.12                 otel_0.2.0                 `\
`#> [25] S4Arrays_1.13.2             tools_4.6.1                `\
`#> [27] parallel_4.6.1              BiocParallel_1.47.0        `\
`#> [29] dplyr_1.2.1                 MultiAssayExperiment_1.39.1`\
`#> [31] curl_8.1.0                  vctrs_0.7.3                `\
`#> [33] R6_2.6.1                    lifecycle_1.0.5            `\
`#> [35] fs_2.1.0                    treeio_1.37.1              `\
`#> [37] pkgconfig_2.0.3             pillar_1.11.1              `\
`#> [39] glue_1.8.1                  Rcpp_1.1.2                 `\
`#> [41] xfun_0.62                   tibble_3.3.1               `\
`#> [43] tidyselect_1.2.1            knitr_1.52                 `\
`#> [45] htmltools_0.5.9             nlme_3.1-171               `\
`#> [47] rmarkdown_2.32              compiler_4.6.1`
