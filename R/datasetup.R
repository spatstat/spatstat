#
#   When the package is installed, this tells us 
#   the directory where the .tab files are stored
#
#   Typically data/murgatroyd.R reads data-raw/murgatroyd.tab
#   and applies special processing
#
spatstat.rawdata.location <- function(...) {
    locn <- system.file("data-raw", package="spatstat")
    if(length(list(...)) != 0) 
      locn <- paste(c(locn, ...), collapse=.Platform$file.sep)
    return(locn)
}
