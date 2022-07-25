#Internal function to parse the attributes/hierarchy list into a data.frame
mgnify_parse_func <- function(json){
    df <- as.data.frame(json$attributes, stringsAsFactors = F)
    df$index_id <- json$attributes$accession
    df
}
