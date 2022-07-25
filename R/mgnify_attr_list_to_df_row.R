# Not exporting this - if people want to they can use the
# rjsonapi functionality. Internally, it takes the "attributes" list
# and converts it into a single row data.frame. For some entries, there is a sublist
# of key/value pairs. metadata_key allows these to be included as columns in the result.
mgnify_attr_list_to_df_row <- function (json, metadata_key=NULL ){

    attrlist=names(json$attributes)

    if (!is.null(metadata_key)){
        baseattrlist <- attrlist[!(attrlist %in% c(metadata_key))]
        metaattrlist <- json$attributes[[metadata_key]]
        metlist <- sapply(metaattrlist, function(x) x$value)
        names(metlist)=sapply(metaattrlist, function(x) x$key)
        df <- as.data.frame(t(unlist(c(json$attributes[baseattrlist], metlist))), stringsAsFactors = F)
    }else{
        df <- as.data.frame(t(unlist(json["attributes"])), stringsAsFactors = F)
    }
    df$accession <- json$id
    df$acc_type <- json$type

    rownames(df) <- df$accession
    df
}
