# Search MGnify database for studies, samples, runs, analyses, biomes, assemblies, and genomes.

Search MGnify database for studies, samples, runs, analyses, biomes,
assemblies, and genomes.

## Usage

``` r
doQuery(x, ...)

# S4 method for class 'MgnifyClient'
doQuery(
  x,
  type = "studies",
  accession = NULL,
  as.df = TRUE,
  max.hits = 200,
  ...
)
```

## Arguments

- x:

  A `MgnifyClient` object.

- ...:

  Remaining parameter key/value pairs may be supplied to filter the
  returned values. Available options differ between `types`. See
  discussion Details section for details.

- type:

  A single character value specifying the type of objects to query. Must
  be one of the following options: `studies`, `samples`, `runs`,
  `analyses`, `biomes`, `assemblies`, `super-studies`,
  `experiment-types`, `pipelines`, `pipeline-tools`, `publications`,
  `genomes`, `genome-search`, `genome-search/gather`,
  `genome-catalogues`, `genomeset`, `cogs`, `kegg-modules`,
  `kegg-classes`, `antismash-geneclusters`, `annotations/go-terms`,
  `annotations/interpro-identifiers`, `annotations/kegg-modules`,
  `annotations/pfam-entries`, `annotations/kegg-orthologs`,
  `annotations/genome-properties`,
  `annotations/antismash-gene-clusters`, `annotations/organisms`, or
  `mydata`. (By default: `type = "studies"`)

- accession:

  A single character value or a vector of character values specifying
  MGnify accession identifiers (of type `type`) or NULL. When NULL, all
  results defined by other parameters are retrieved. (By default:
  `accession = NULL`)

- as.df:

  A single boolean value specifying whether to return the results as a
  data.frame or leave as a nested list. In most cases, `as.df = TRUE`
  will make the most sense. (By default: `as.df = TRUE`)

- max.hits:

  A single integer value specifying the maximum number of results to
  return or FALSE. The actual number of results will actually be higher
  than `max.hits`, as clipping only occurs on pagination page
  boundaries. To disable the limit, set `max.hits = NULL`. (By default:
  `max.hits = 200`)

## Value

A nested list or data.frame containing the results of the query.

## Details

`doQuery` is a flexible query function, harnessing the "full" power of
the JSONAPI MGnify search filters. Search results may be filtered by
metadata value, associated study/sample/analyse etc.

See [Api browser](https://www.ebi.ac.uk/metagenomics/api/v1/) for
information on MGnify database filters. You can find help on customizing
queries from
[here](https://emg-docs.readthedocs.io/en/latest/api.html#customising-queries).

For example the following filters are available:

- **studies**: accession, biome_name, lineage, centre_name, include

- **samples**: accession, experiment_type, biome_name, lineage,
  geo_loc_name, latitude_gte, latitude_lte, longitude_gte,
  longitude_lte, species, instrument_model, instrument_platform,
  metadata_key, metadata_value_gte, metadata_value_lte, metadata_value,
  environment_material, environment_feature, study_accession, include

- **runs**: accession, experiment_type, biome_name, lineage, species,
  instrument_platform, instrument_model, metdata_key,
  metadata_value_gte, metadata_value_lte, metadata_value,
  sample_accession, study_accession, include

- **analyses**: biome_name, lineage, experiment_type, species,
  sample_accession, pipeline_version

- **biomes**: depth_gte, depth_lte

- **assemblies**: depth_gte, depth_lte

Unfortunately it appears that in some cases, some of these filters don't
work as expected, so it is important to check the results returned match
up with what's expected. Even more unfortunately if there's an error in
the parameter specification, the query will run as if no filter
parameters were present at all. Thus the result will appear
superficially correct but will infact correspond to something completely
different. This behaviour will hopefully be fixed in future incarnations
of the MGnifyR or JSONAPI, but for now users should double check
returned values.

It is currently not possible to combine queries of the same type in a
single call (for example to search for samples *between* latitude).
However, it is possible to run multiple queries and combine the results
using set operations in R to get the desired behaviour.

## Examples

``` r
mg <- MgnifyClient(useCache = FALSE)

# Get a list of studies from the Agricultural Wastewater :
agwaste_studies <- doQuery(
    mg, "studies", biome_name="Agricultural wastewater"
    )

if (FALSE) { # \dontrun{
# Get all samples from a particular study
samps <- doQuery(mg, "samples", accession="MGYS00004521")

# Search polar samples
samps_np <- doQuery(mg, "samples", latitude_gte=66, max.hits=10)
samps_sp <- doQuery(mg, "samples", latitude_lte=-66, max.hits=10)

# Search studies that have studied drinking water
tbl <- doQuery(
    mg,
    type = "studies",
    biome_name = "root:Environmental:Aquatic:Freshwater:Drinking water",
    max.hits = 10)
} # }
```
