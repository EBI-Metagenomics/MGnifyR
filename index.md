# Metagenomics workshop MGnifyR session details.

## Introduction

The exercises here don't require a particularly high level of proficiency in R, but you will need to understand basic R concepts such as calling functions with parameters, and basic data.frame manipulation to complete the steps. It might help to review a tutorial like http://www.r-tutor.com/r-introduction/data-frame before you start. That said, the intention is not to make all students experts with R and `MGnifyR` by the end of the session, not least because the package still in development and is liable to change. Instead, this is an opportunity for non-developer users to have a go with `MGnifyR` and provide feedback about what works and what doesn't. What features do you like? What would you like to see in the future?.

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
 


## Exercises:

### Basics

 - Start up an R session in whatever tool you feel most comfortable with. RStudio, R-gui and command line are all available on the virtual machine image.
 
 - Using Devtools and `install_github("beadyallen/MGnifyR")` command, update the version MGnifyR package installed on your VM. Installation and updating use exactly the same commands. Since the image was created, numerous bugs have been identified so it's worth doing the update before we begin. No other packages will require updating, when asked you may safely decline updating the other packages. 
 
 - After installation, load the library into your session: `library(MGnfifyR)`
 
 - Create your session `mgnify_client` object, you'll need it for the other exercises.
    Making sure you've enabled the local result cache, where will it be stored? Open up a file browser and check to see the directory has been created properly. 
    

### Searching MGnify

 - Using the `mgnify_query` command, search for all lentic samples - just as you did this morning with the web interface. You'll need the `biome_name` parameter in your search, and the full colon seperated biome path. 
    - How many samples can you find? Why is the count different from this morning? Can you restrict the results to just those from this morning?
    - Subset the `data.frame` to keep only samples with depths between 25m and 50m. How many are there? What extra step do you need to do to make this work?   
    
    
### Taxa abundances with `phyloseq`

 - Retrieve the complete metadata for studies MGYS00002662 (Anaerobic digestion from York University), and MGYS00001910 (an activated sludge retention time survey from EAWAG). 
    - Which experiment types are included in these studies.
    - How many samples are present from each study?
    
 - Using `mgnify_get_analyses_phyloseq`, build the `phyloseq` object for the amplicon experiment accessions only.
    - How many unique taxa are present across all samples in both studies? 
    - Plot normalized taxa abundances across all samples using the `plot_bar` function in `phyloseq`. Using the vignette as a guide, 
    - Which 
    
  MGYS00002662 - York uni anaerobic digestion amplicons
  MGYS00001910 - EAWAG Activated sludge
  
 
 
 
### Human skin functional results
 - Earlier you looked at a single analysis accession from a larger metagenomics study. With the power of `MGnifyR`, retrieve the full sample, study and analysis metadata `data.frame` for study **MGYS00003598**. Remember that result retrieval is based on analyses, so the first step is to get the set of associated analysis accessions.
    - According to MGnifyR, how many samples have been included in this study? 
    - Are you sure? Which pipeline versions are present?
    
 - Using only the results for the latest pipeline, retrieve the go-slim annotation for 
     - For sample 
     
     
 - Get downloads for assembly MGYA00510849:
    -Ret
 

 assembly accession is "MGYA00510849"
 
 ## Exercise hints
 
 ### Searching:
   - The `biome_name` to search for is `root:Environmental:Aquatic:Lentic`.
   - You might use the R `table` command to count up unique occurrences in a column or vector.
   - **ALL** metadata columns are "characters" by default. If you was to do numeric comparisons, you'll have to convert to numbers. Use `as.numeric`.
   
   
 ### Phyloseq data
 
 
 ### Functional results

