Load every single function in R/Mgnify.R

    devtools::load_all()

    ## ℹ Loading MGnifyR

    library(MGnifyR)

    # Set up the MGnify client instance

    mgclnt <- mgnify_client(usecache = T, cache_dir = '/tmp/MGnify_cache')

    # Retrieve the list of analyses associated with a study

    accession_list <- mgnify_analyses_from_studies(mgclnt, "MGYS00005058", usecache = T)

    ##   |                                                                                                                                                                             |                                                                                                                                                                     |   0%  |                                                                                                                                                                             |=====================================================================================================================================================================| 100%

    # Download all associated study/sample and analysis metadata

    meta_dataframe <- mgnify_get_analyses_metadata(mgclnt, accession_list, usecache = T )

    ##   |                                                                                                                                                                             |                                                                                                                                                                     |   0%  |                                                                                                                                                                             |=======                                                                                                                                                              |   4%  |                                                                                                                                                                             |==============                                                                                                                                                       |   8%  |                                                                                                                                                                             |=====================                                                                                                                                                |  12%  |                                                                                                                                                                             |============================                                                                                                                                         |  17%  |                                                                                                                                                                             |==================================                                                                                                                                   |  21%  |                                                                                                                                                                             |=========================================                                                                                                                            |  25%  |                                                                                                                                                                             |================================================                                                                                                                     |  29%  |                                                                                                                                                                             |=======================================================                                                                                                              |  33%  |                                                                                                                                                                             |==============================================================                                                                                                       |  38%  |                                                                                                                                                                             |=====================================================================                                                                                                |  42%  |                                                                                                                                                                             |============================================================================                                                                                         |  46%  |                                                                                                                                                                             |==================================================================================                                                                                   |  50%  |                                                                                                                                                                             |=========================================================================================                                                                            |  54%  |                                                                                                                                                                             |================================================================================================                                                                     |  58%  |                                                                                                                                                                             |=======================================================================================================                                                              |  62%  |                                                                                                                                                                             |==============================================================================================================                                                       |  67%  |                                                                                                                                                                             |=====================================================================================================================                                                |  71%  |                                                                                                                                                                             |============================================================================================================================                                         |  75%  |                                                                                                                                                                             |===================================================================================================================================                                  |  79%  |                                                                                                                                                                             |==========================================================================================================================================                           |  83%  |                                                                                                                                                                             |================================================================================================================================================                     |  88%  |                                                                                                                                                                             |=======================================================================================================================================================              |  92%  |                                                                                                                                                                             |==============================================================================================================================================================       |  96%  |                                                                                                                                                                             |=====================================================================================================================================================================| 100%

    # Convert analyses outputs to a single `TreeSummarizedExperiment` object
    tse <- mgnify_get_analyses_treese(mgclnt, meta_dataframe$analysis_accession, usecache = T)

    ##   |                                                                                                                                                                             |                                                                                                                                                                     |   0%  |                                                                                                                                                                             |=======                                                                                                                                                              |   4%  |                                                                                                                                                                             |==============                                                                                                                                                       |   8%  |                                                                                                                                                                             |=====================                                                                                                                                                |  12%  |                                                                                                                                                                             |============================                                                                                                                                         |  17%  |                                                                                                                                                                             |==================================                                                                                                                                   |  21%  |                                                                                                                                                                             |=========================================                                                                                                                            |  25%  |                                                                                                                                                                             |================================================                                                                                                                     |  29%  |                                                                                                                                                                             |=======================================================                                                                                                              |  33%  |                                                                                                                                                                             |==============================================================                                                                                                       |  38%  |                                                                                                                                                                             |=====================================================================                                                                                                |  42%  |                                                                                                                                                                             |============================================================================                                                                                         |  46%  |                                                                                                                                                                             |==================================================================================                                                                                   |  50%  |                                                                                                                                                                             |=========================================================================================                                                                            |  54%  |                                                                                                                                                                             |================================================================================================                                                                     |  58%  |                                                                                                                                                                             |=======================================================================================================                                                              |  62%  |                                                                                                                                                                             |==============================================================================================================                                                       |  67%  |                                                                                                                                                                             |=====================================================================================================================                                                |  71%  |                                                                                                                                                                             |============================================================================================================================                                         |  75%  |                                                                                                                                                                             |===================================================================================================================================                                  |  79%  |                                                                                                                                                                             |==========================================================================================================================================                           |  83%  |                                                                                                                                                                             |================================================================================================================================================                     |  88%  |                                                                                                                                                                             |=======================================================================================================================================================              |  92%  |                                                                                                                                                                             |==============================================================================================================================================================       |  96%  |                                                                                                                                                                             |=====================================================================================================================================================================| 100%

    # retrievelist requires us to specify which functional estimations we

    # wish to examine, and may be one or more of go-slim, go-terms,

    # interpro-identifiers or antismash-gene-clusters. retrievelist may

    # also be all, in which case all available results are retrieved,

    # along with all taxonomic assignments.

    # Retrieve assignment counts for these analyses as a data list

    #ip_df <- mgnify_get_analyses_results(mgclnt, meta_dataframe$analysis_accession, retrievelist = c("interpro-identifiers"), usecache = T)

    dli <- mgnify_get_analyses_results(mgclnt, meta_dataframe$analysis_accession, retrievelist = "all", usecache = T)

    ##   |                                                                                                                                                                             |                                                                                                                                                                     |   0%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |=======                                                                                                                                                              |   4%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |==============                                                                                                                                                       |   8%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |=====================                                                                                                                                                |  12%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |============================                                                                                                                                         |  17%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |==================================                                                                                                                                   |  21%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |=========================================                                                                                                                            |  25%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |================================================                                                                                                                     |  29%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |=======================================================                                                                                                              |  33%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |==============================================================                                                                                                       |  38%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |=====================================================================                                                                                                |  42%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |============================================================================                                                                                         |  46%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |==================================================================================                                                                                   |  50%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |=========================================================================================                                                                            |  54%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |================================================================================================                                                                     |  58%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |=======================================================================================================                                                              |  62%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |==============================================================================================================                                                       |  67%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |=====================================================================================================================                                                |  71%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |============================================================================================================================                                         |  75%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |===================================================================================================================================                                  |  79%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |==========================================================================================================================================                           |  83%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |================================================================================================================================================                     |  88%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |=======================================================================================================================================================              |  92%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |==============================================================================================================================================================       |  96%

    ## Warning: Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.
    ## Unknown or uninitialised column: `index_id`.

    ##   |                                                                                                                                                                             |=====================================================================================================================================================================| 100%

    # Seems that taxonomy-ssu corresponds to the phyloseq object

    # all(abundances(psobj) == dli[["taxonomy-ssu"]][, 9:32])

    # TODO

    # 1. mgnify_get_analyses_phyloseq -> mgnify_get_analyses_treese (new function)

    # 2. mgnify_get_analyses_results -> mgnify_get_analyses_results_mae (new function; to download & arrange selected tables in a MultiAssayExperiment)

    # 3. some use examples in MGnifyR vignette

    # 4. some use examples in OMA

    # 5. MGnifyR to CRAN (joint project?)
