#Internal function to parse the attributes/hierarchy list into a data.frame
mgnify_parse_tax <- function(json){
  df <- as.data.frame(c(json$attributes["count"], unlist(json$attributes$hierarchy)), stringsAsFactors = F)
  df$index_id <- json$attributes$lineage
  df
}
