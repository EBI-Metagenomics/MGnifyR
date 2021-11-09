# Metagenomics workshop MGnifyR practical session.

# Introduction

Most of the exercises here don't require a particularly high level of proficiency in R, but you will need to understand basic R concepts such as calling functions with parameters, and basic data.frame manipulation to complete the steps. It might help to review a tutorial like http://www.r-tutor.com/r-introduction/data-frame before you start. That said, the intention is not to make all students experts with R and `MGnifyR` by the end of the session, not least because the package still in development and is liable to change. Instead, this is an opportunity for non-developer users to have a go with `MGnifyR` and provide feedback about what works and what doesn't. What features do you like? What would you like to see in the future?.

Detailed help for each function is available in R using the standard "?function_name" command (i.e. typing `?mgnify_query` will bring up built-in help for the `mgnify_query` command), and a vignette is available containing a reasonably verbose overview of the main functionality. This can be read either within R with the vignette("MGnifyR") command, or from the development repository at https://htmlpreview.github.io/?https://github.com/beadyallen/MGnifyR/blob/master/doc/MGnifyR.html

### Command cheat sheet
The following list of key functions should give a starting point for finding relevent documentation.

 - `mgnify_client()` : Create the client object required for all other functions.
 - `mgnify_query()` : Search the whole MGnify database
 - `mgnify_analyses_from_xxx()` : Convert `xxx` accessions to `analyses` accessions. `xxx` is either `samples` or `studies`
 - `mgnify_get_analyses_metadata()` : Retrieve all study, sample and analysis metadata for given analyses.
 - `mgnify_get_analyses_phyloseq()` : Convert abundance, taxonomic, and sample metadata into a single `phyloseq` object.
 - `mgnify_get_analyses_results()` : Get functional annotation results for a set of analyses.
 - `mgnify_download()` : Download raw results files from MGnify
 - `mgnify_retrieve_json()` : Low level API access helper function.
 


# Exercises:

### Basics
 - Start up an R session in whatever tool you feel most comfortable with. Both RStudio and command line are available on the virtual machine image.
 
 - Using Devtools and `install_github("beadyallen/MGnifyR")` command, update the version MGnifyR package installed on your VM. Installation and updating use exactly the same commands. Since the image was created, numerous bugs have been identified so it's worth doing the update before we begin. No other packages will require updating, when asked you may safely decline updating the other packages. 
 
 - After installation, load the library into your session: `library(MGnfifyR)`
 
 - Create your session `mgnify_client` object, you'll need it for the other exercises.
     - Try enabling the caching features. Where will the files be stored? Open up a file browser and check to see the directory has been created properly.
    
### Searching MGnify
 - Using the `mgnify_query` command, search for all lentic samples - just as you did this morning with the web interface. You'll need the `biome_name` parameter in your search, and the full colon seperated biome path. 
    - How many samples can you find? Why is the count different from this morning? Can you restrict the results to just those from this morning?
    - Subset the `data.frame` to keep only samples with depths between 25m and 50m. How many are there? What extra step do you need to do to make this work?   
    
### Taxa abundances with `phyloseq`
 - Retrieve the complete metadata for studies **MGYS00002662** (Anaerobic digestion from York University), and **MGYS00001910** (an activated-sludge sludge retention time survey from EAWAG). 
    - Which experiment types are included in these studies.
    - How many analyses are present in each study?
    
 - Using `mgnify_get_analyses_phyloseq`, build the `phyloseq` object *for the amplicon* experiment accessions only.
    - Optionally, try filtering the analyses, removing those with low total abundances. A `hist` plot can be useful to determine a good threshold.
    - How many unique taxa are present across all samples in both studies? 
    - Plot normalized taxa abundances across all samples using `plot_bar`. 
    - Using `phyloseq`'s `estimate_richness`, which sample has the highest taxonomic diversity (by observed species)? Can you produce a box plot of species richness grouped by source biome? Follow the example in the vignette for help with the plotting.
    - If feeling adventurous, try producing an ordination plot as well displaying the taxonomic separation between environment types. `metaMDS` in the `vegan` package will work but will require a reasonable amount of data munging. An example is shown in the `MGnifyR` vignette, although not with `phyloseq`.
   
### Human skin functional results
 - Earlier today you looked at a single analysis accession from a larger metagenomics study. With the power of `MGnifyR`, retrieve the full sample, study and analysis metadata `data.frame` for study **MGYS00003598**. Remember that result retrieval is based on analyses, so the first step is to get the set of associated analysis accessions.
    - According to MGnifyR, how many `samples` have been included in this study? 
    - Are you sure? Which pipeline versions are present?
    - Which analysis contains the greatest number of pCDS sequences? What about the analysis with the most Interpro annotations?
    
 - Using only the results for the most up-to-date pipeline version, retrieve the corresponding go-slim annotations.
     - Which GO term is most abundant across all samples (by mean abundance per sample)?
     
 - Get the Interpro identifiers for the single sample identified above. 
     - For this sample, which identifier is most abundant?
     - Does this agree with the GO terms result?
 
### Downloading "raw" files.
   - Sticking with the functional results from the previous step, determine the URL locations of *all* downloadable files for the analyses you found with pipeline version 5.0. 
      - How many files are available for download across the version 5.0 analyses?
      - How many different types of output are there?
      - Try downloading a file of your choice from those available. Try navigating to it afterwards in the file explorer. 
      
   - For a final, advanced exercise, download both the protein sequence FASTA file and Interpro annotations for analysis MGYA00510855. Find the complete set of unique protein sequences annotated with interpro id IPR036388 (the most abundant protein). The last section of the vignette will probably prove useful to perform this task. 
      - How many unique sequences are there? Does it match what we found earlier from the Interpro results?
   
   
   
# Exercise hints - have a look if you're getting stuck.
 
   - Use `colnames` to list data.frame columns. MGnifyR results can be quite unwieldy to view all at once, so subsetting only those columns of interest can be useful.
 
### Searching:
   - The `biome_name` to search for is `root:Environmental:Aquatic:Lentic`.
   - You might use the R `table` command to count up unique occurrences in a column or vector.
   - **ALL** metadata columns are "characters" by default. If you was to do numeric comparisons, you'll have to convert to numbers. Use `as.numeric`.
   
### Phyloseq data
   - The terminology can get a little mixed up. `samples` in phyloseq really correspond to `analyses` in MGnifyR. 
   - Use `mgnify_analyses_from_studies` to get an accession list, then `mgnify_get_analysis_metadata` to retrieve the metadata.
   - This section might be a bit daunting, and will be unfamiliar if you've not got much R experience. Don't worry - retrieving the phyloseq object is the main objective. Once you have it in this format, there's loads of tutorials and walkthroughs available on the internet. The [phyloseq](https://joey711.github.io/phyloseq/)  website is a great place to start.
   - Just displaying the downloaded object should get you the answer to the first question.
   - There's many approaches to cross sample normalization, but although it's frowned upon by some, `phyloseq`'s `rarefy_even_depth` will work fine in this case - provided low abundance samples are removed.
   - Useful `phyloseq` functions include `sample_sums` (total abundance per-sample) and `subset_samples` (remove samples from the `phyloseq` object based on a logical vector). The OTU abundance table is accessed with `otu_table()`,  and the taxonomy assignment table can be extracted from the phyloseq object with `tax_table`.
 
### Functional results
   - Convert to an `analysis` accession list with `mgnify_analyses_from_studies`. then use `mgnify_get_analyses_metadata`.
   - The first few questions really only concern the metadata. 
   - Don't assume all the rows are unique - this will only be true for the `analysis` accession column.
   - If you have trouble with finding which particular entry is the most abundant, don't worry. Just pick a valid id from the results you've got and try the next step. Mastering R data munging is *not* the point of the exercise.
   - Use `mgnify_get_analyses_results` with a suitable setting for the `retrievelist` parameter. Use the online help to check what you need.
   - Function results are always returned as a `list` of `data.frames`.
   - Try `rowMeans` - making sure to drop the non-numeric columns (indexing a dataframe with a -ve value removes rather than selects a row/column).
 
