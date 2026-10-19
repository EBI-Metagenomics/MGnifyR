# Metagenomics bioinformatics at MGnify

## Introduction

In this notebook we aim to demonstrate how the `MGnifyR` tool can be
used to fetch data of a MGnify microbiome data resource. We then
showcase how to analyze the datausing advanced microbiome data science
tools, including estimating alpha and beta diversity, as well as
performing differential abundance analysis.

[`MGnifyR`](https://www.bioconductor.org/packages/release/bioc/html/MGnifyR.html)
is an R/Bioconductor package that provides a set of tools for easily
accessing and processing MGnify data in R, making queries to MGnify
databases through the [MGnify
API](https://www.ebi.ac.uk/metagenomics/api/v1/).

The benefit of `MGnifyR` is that it streamlines data access, allowing
users to fetch data either in its “raw” format or directly as a
[`TreeSummarizedExperiment`
(`TreeSE`)](https://microbiome.github.io/OMA/docs/devel/pages/containers.html)
object. This enables seamless integration with custom workflows for
analysis.

Utilizing `TreeSE` provides access to a wide range of tools within
Bioconductor’s `SummarizedExperiment` (`SE`) ecosystem. It also
integrates with the [`mia` package](https://microbiome.github.io/mia/),
which offers microbiome-specific methods within the `SE` framework.

For more information on microbiome data science in Bioconductor, refer
to [Orchestrating Microbiome Analysis (OMA) online
book](https://microbiome.github.io/OMA/docs/devel/).

## Load packages

\
`# List of packages that we need`\
`packages`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``"ANCOMBC"``, ``"MGnifyR"``, ``"mia"``,  ``"miaViz"``, ``"scater"``)`\
\
`# Get packages that are already installed`\
`packages_already_installed`` ``<-`` ``packages``[`` ``packages`` `[`%in%`](https://rdrr.io/pkg/BiocGenerics/man/match.html)` `[`installed.packages`](https://rdrr.io/r/utils/installed.packages.html)`(``)`` ``]`\
`# Get packages that need to be installed`\
`packages_need_to_install`` ``<-`` `[`setdiff`](https://generics.r-lib.org/reference/setops.html)`(`` ``packages``, ``packages_already_installed`` ``)`\
`# Loads BiocManager into the session. Install it if it is not already installed.`\
`if``(`` ``!`[`require`](https://rdrr.io/r/base/library.html)`(`[`"BiocManager"`](https://bioconductor.github.io/BiocManager/)`)`` ``)``{`\
`    `[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"BiocManager"``)`\
`    `[`library`](https://rdrr.io/r/base/library.html)`(`[`"BiocManager"`](https://bioconductor.github.io/BiocManager/)`)`\
`}`\
`# If there are packages that need to be installed, installs them with BiocManager`\
`# Updates old packages.`\
`if``(`` `[`length`](https://rdrr.io/r/base/length.html)`(``packages_need_to_install``)`` ``>`` ``0`` ``)`` ``{`\
`   `[`install`](https://bioconductor.github.io/BiocManager/reference/install.html)`(``packages_need_to_install``, ask ``=`` ``FALSE``)`\
`}`\
\
`# Load all packages into session. Stop if there are packages that were not`\
`# successfully loaded`\
`pkgs_not_loaded`` ``<-`` ``!`[`sapply`](https://rdrr.io/pkg/BiocGenerics/man/lapply.html)`(``packages``, ``require``, character.only ``=`` ``TRUE``)`\
`pkgs_not_loaded`` ``<-`` `[`names`](https://rdrr.io/r/base/names.html)`(``pkgs_not_loaded``)``[`` ``pkgs_not_loaded`` ``]`\
`if``(`` `[`length`](https://rdrr.io/r/base/length.html)`(``pkgs_not_loaded``)`` ``>`` ``0`` ``)``{`\
`    `[`stop`](https://rdrr.io/r/base/stop.html)`(`\
`        ``"Error in loading the following packages into the session: '"``,`\
`        `[`paste0`](https://rdrr.io/r/base/paste.html)`(``pkgs_not_loaded``, collapse ``=`` ``"', '"``)``, ``"'"``)`\
`}`

## Data import

To interact with the MGnify database, we need to create a `MgnifyClient`
object. This object allows us to store options for data fetching. For
instance, we can configure it to use a cache for improved efficiency.

\
`# Create the MgnifyClient object with caching enabled`\
`mg`` ``<-`` `[`MgnifyClient`](../reference/MgnifyClient.md)`(`\
`  useCache ``=`` ``TRUE``,`\
`  cacheDir ``=`` ``"/home/training"`` ``# Set this to your desired cache directory`\
`  ``)`

In this workflow, we will fetch taxonomy annotations and metadata from
the study
[“MGYS00005154”](https://www.ebi.ac.uk/metagenomics/studies/MGYS00005154).
The dataset focuses on the human gut microbiome, analyzed across
different geographic regions.

We can now search for all analyses associated with the certain study.
The analysis refers to metagenomic runs performed to samples. Each
sample can have multiple runs made, which is why we work with analyses
and not with samples; analysis identifier points to a single entity.

\
`study_id`` ``<-`` ``"MGYS00005154"`\
`analysis_id`` ``<-`` `[`searchAnalysis`](../reference/searchAnalysis.md)`(``mg``, ``"studies"``, ``study_id``)`

Now we are ready to load the metadata on the analyses to get idea on
what kind of data we are dealing with.

There are currently (17 Sep 2024) almost 1,000 analyses available.
Downloading whole dataset will take some time, which is why we use
memory cache.

\
`metadata`` ``<-`` `[`getMetadata`](../reference/getMetadata.md)`(``mg``, accession ``=`` ``analysis_id``)`

We can see that there are analyses that are performed with different
pipelines. Let’s take only those analyses that are generated with the
pipeline version 5.0.

\
`metadata`` ``<-`` ``metadata``[``metadata``[[``"analysis_pipeline-version"``]``]`` ``==`` ``"5.0"``, ``]`

We have now analyses that each point to unique sample. The final step is
to fetch abundance tables in `TreeSummarizedExperiment` (`TreeSE`)
format.

\
`tse`` ``<-`` `[`getResult`](../reference/getResult.md)`(`\
`    ``mg``,`\
`    accession ``=`` ``metadata``[[``"analysis_accession"``]``]``,`\
`    get.func ``=`` ``FALSE`\
`    ``)`\
`tse`

The fetched data is `TreeSE` object, including taxonomy annotations. See
[OMA online
book](https://microbiome.github.io/OMA/docs/devel/pages/containers.html)
on how to handle the data in this format.

## Preprocessing

Below, we agglomerate the data to the Order level, meaning we summarize
the abundances at this specific taxonomic rank. The OMA provides a
detailed
[chapter](https://microbiome.github.io/OMA/docs/devel/pages/agglomeration.html)
explaining agglomeration in more depth.

\
`tse_order`` ``<-`` ``agglomerateByRank``(``tse``, rank ``=`` ``"Order"``)`

Because of the unique properties of microbiome data, we have to apply
transformations. Here, we perform relative transformation. You can find
more information on transformations from
[OMA](https://microbiome.github.io/OMA/docs/devel/pages/transformation.html).

\
`# Transform the main TreeSE`\
`tse`` ``<-`` ``transformAssay``(``tse``, method ``=`` ``"relabundance"``)`\
`# Transform the agglomerated TreeSE`\
`tse_order`` ``<-`` ``transformAssay``(``tse_order``, method ``=`` ``"relabundance"``)`

## Alpha diversity

Alpha diversity measures community diversity within a sample. Learn more
on community diversity from
[here](https://microbiome.github.io/OMA/docs/devel/pages/alpha_diversity.html).

\
`# Calculate alpha diversity`\
`tse`` ``<-`` ``addAlpha``(``tse``)`\
\
`# Create a plot`\
`p`` ``<-`` ``plotColData``(`\
`  ``tse``,`\
`  y ``=`` ``"shannon_diversity"``,`\
`  x ``=`` ``"sample_geographic.location..country.and.or.sea.region."``,`\
`  show_boxplot ``=`` ``TRUE`\
`  ``)`\
`p`

We can test whether the diversity differences are statistically
significant. We utilize Mann Whithney U test (or Wilcoxon test).

\
[`pairwise.wilcox.test`](https://rdrr.io/r/stats/pairwise.wilcox.test.html)`(`\
`    ``tse``[[``"shannon_diversity"``]``]``,`\
`    ``tse``[[``"sample_geographic.location..country.and.or.sea.region."``]``]``,`\
`    p.adjust.method ``=`` ``"fdr"`\
`    ``)`

To add p-values to the plot, see
[OMA](https://microbiome.github.io/OMA/docs/devel/pages/alpha_diversity.html#visualizing-significance-in-group-wise-comparisons).

## Beta diversity

We can assess the differences in microbial compositions between samples,
aiming to identify patterns in the data that are associated with
covariates.

To achieve this, we perform Principal Coordinate Analysis (PCoA) using
Bray-Curtis dissimilarity.

\
`# Perform PCoA`\
`tse`` ``<-`` ``runMDS``(`\
`  ``tse``,`\
`  FUN ``=`` ``getDissimilarity``,`\
`  method ``=`` ``"bray"``,`\
`  assay.type ``=`` ``"relabundance"`\
`)`\
`# Visualize PCoA`\
`p`` ``<-`` ``plotReducedDim``(`\
`  ``tse``,`\
`  dimred ``=`` ``"MDS"``,`\
`  colour_by ``=`` ``"sample_geographic.location..country.and.or.sea.region."`\
`  ``)`\
`p`

See [community similarity
chapter](https://microbiome.github.io/OMA/docs/devel/pages/beta_diversity.html)
from OMA for more information.

## Differential abundance analysis (DAA)

In DAA, we analyze whether abundances of certain features vary between
study groups. Again, OMA has a dedicated chapter also on this
[topic](https://microbiome.github.io/OMA/docs/devel/pages/differential_abundance.html).

\
`# Perform DAA`\
`res`` ``<-`` ``ancombc2``(`\
`    data ``=`` ``tse_order``,`\
`    assay.type ``=`` ``"counts"``,`\
`    fix_formula ``=`` ``"sample_geographic.location..country.and.or.sea.region."``,`\
`    p_adj_method ``=`` ``"fdr"``,`\
`    ``)`

Next we visualize features that have the lowest p-values.

\
`# Get the most significant features`\
`n_top`` ``<-`` ``4`\
`res_tab`` ``<-`` ``res``[[``"res"``]``]`\
`res_tab`` ``<-`` ``res_tab``[`[`order`](https://rdrr.io/pkg/BiocGenerics/man/order.html)`(``res_tab``[[``"q_(Intercept)"``]``]``)``, ``]`\
`top_feat`` ``<-`` ``res_tab``[`[`seq_len`](https://rdrr.io/r/base/seq.html)`(``n_top``)``, ``"taxon"``]`\
\
`# Create a plot`\
`p`` ``<-`` ``plotExpression``(`\
`  ``tse_order``,`\
`  features ``=`` ``top_feat``,`\
`  assay.type ``=`` ``"relabundance"``,`\
`  x ``=`` ``"sample_geographic.location..country.and.or.sea.region."``,`\
`  show_boxplot ``=`` ``TRUE``, show_violin ``=`` ``FALSE``, point_shape ``=`` ``NA`\
`  ``)`` ``+`\
`  ``scale_y_log10``(``)`\
`p`

## Session info

\
[`sessionInfo`](https://rdrr.io/r/utils/sessionInfo.html)`(``)`
