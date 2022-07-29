#Retrieves combined study/sample/analysis metadata - not exported
mgnify_get_single_analysis_metadata <- function(client=NULL, accession, usecache=T, maxhits=-1){

    dat <- mgnify_retrieve_json(client, paste("analyses", accession, sep="/"), usecache = usecache, maxhits = maxhits)
    #There ~should~ be just a single result
    top_data <- dat[[1]]
    analysis_df <- mgnify_attr_list_to_df_row(top_data, metadata_key = "analysis-summary")

    #cat(str(analysis_metadata))
    #Build up the metadata dataframe from the analyses_metadata_headers vector:
    sample_met <- mgnify_retrieve_json(client, complete_url = top_data$relationships$sample$links$related, usecache = usecache)
    study_met <- mgnify_retrieve_json(client, complete_url = top_data$relationships$study$links$related, usecache = usecache)

    sample_df <- mgnify_attr_list_to_df_row(sample_met[[1]], metadata_key = "sample-metadata")
    #It turns out that a sample might not be part of a study - if it's been harvested...
    #So tryCatch it and return an empy df row if things go south.
    study_df <- tryCatch(mgnify_attr_list_to_df_row(study_met[[1]]), error=function(X) {
        warning(paste("Failed to find study metadata for ", accession, sep=""))
        data.frame(accession=NA)
    })

    colnames(sample_df) <- paste("sample",colnames(sample_df), sep="_")
    colnames(study_df) <- paste("study",colnames(study_df), sep="_")
    colnames(analysis_df) <- paste("analysis",colnames(analysis_df), sep="_")

    rownames(sample_df) <- rownames(analysis_df)
    rownames(study_df) <- rownames(analysis_df)
    full_df <- cbind(analysis_df, study_df, sample_df)

    #extras - include some more metadata from various places
    #assembly accesion
    if("id" %in% names(top_data$relationships$assembly$data)){
        full_df$assembly_accession <- top_data$relationships$assembly$data$id
    }
    #run accession
    if("id" %in% names(top_data$relationships$run$data)){
        full_df$run_accession <- top_data$relationships$run$data$id
    }

    #biom (from the sample metadata)
    tryCatch({
        full_df$biome_string <- sample_met[[1]]$relationships$biome$data$id
    },
        error <- function(x) warning("Error finding biome entry")
    )

    full_df
}

