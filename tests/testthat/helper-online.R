# Helpers for tests that use the servers: MGnify API, the FTP server of
# MGnify, where the files are, and BioSamples. Problems of the network or the
# servers, e.g., a server not responding, are not problems of the package.
# Such tests are skipped instead of failed, so that automated testing does
# not give false alarms.

# The servers, each with a small request that should always work
.TEST_SERVERS <- c(
    "MGnify API" =
        "https://www.ebi.ac.uk/metagenomics/api/v2/studies/?page_size=1",
    "MGnify FTP" = "https://ftp.ebi.ac.uk/pub/databases/metagenomics/",
    "BioSamples" = paste0(
        "https://www.ebi.ac.uk/biosamples/v2/samples/bulk-fetch?",
        "accessions=SAMN02595037"))

# Get the servers that do not respond. A server is tried 3 times, because the
# FTP server, for example, can refuse single connections while it works
# otherwise.
servers_down <- function(){
    responds <- vapply(.TEST_SERVERS, function(url){
        for( attempt in seq_len(3L) ){
            handle <- curl::new_handle(timeout = 30)
            res <- tryCatch(
                curl::curl_fetch_memory(url, handle = handle),
                error = function(e) NULL)
            if( !is.null(res) && res[["status_code"]] < 500 ){
                return(TRUE)
            }
            Sys.sleep(attempt)
        }
        return(FALSE)
    }, logical(1))
    down <- names(.TEST_SERVERS)[!responds]
    return(down)
}

# Run a test that uses the servers:
# 1. The test is skipped if any server is down when the test starts.
# 2. If the test fails, the servers are checked again, and the test is
#    skipped if any of them is down, because the failure was then likely
#    caused by the server.
# 3. If the servers respond, the test is run once more, because a server can
#    fail only for a moment, e.g., one download fails. The data that was
#    already fetched is in the cache, so the second run is fast. A failure in
#    the second run is a failure of the test.
test_online <- function(desc, code){
    code <- substitute(code)
    env <- parent.frame()
    test_that(desc, {
        down <- servers_down()
        if( length(down) > 0L ){
            skip(paste("Not responding:", paste(down, collapse = ", ")))
        }
        # Failed expectations are also errors, so they are caught here
        first <- tryCatch({
            eval(code, new.env(parent = env))
            NULL
        }, error = function(e) e)
        if( !is.null(first) ){
            down <- servers_down()
            if( length(down) > 0L ){
                skip(paste(
                    "Stopped responding during the test:",
                    paste(down, collapse = ", "), "-",
                    conditionMessage(first)))
            }
            eval(code, new.env(parent = env))
        }
    })
}
