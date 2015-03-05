##  copyExampleFiles.R
##  $Revision: 1.7 $ $Date: 2015/03/05 11:42:01 $

copyExampleFiles <- function(which, folder=getwd()) {
  choices <- dir(system.file("rawdata", package="spatstat"))
  if(missing(which) || is.null(which)) {
    message(paste("Choices are: which=", commasep(sQuote(choices), " or ")))
    return(invisible(NULL))
  }
  whichdata <- match.arg(which, choices)
  sourcefolder <- system.file("rawdata", whichdata, package="spatstat")
  sourcefiles <- dir(sourcefolder)
  if(length(sourcefiles) == 0)
    stop("No files available")
  # set directory
  oldfolder <- getwd()
  setwd(folder)
  on.exit(setwd(oldfolder))
  # 
  for(fn in sourcefiles) {
    frompath <- file.path(sourcefolder, fn)
    file.copy(from = frompath, to = fn, overwrite=TRUE)
  }
  splat("Copied",
        ngettext(length(sourcefiles), "file", "files"),
        commasep(dQuote(sourcefiles)), "to",
        if(identical(folder, oldfolder)) "current folder" else dQuote(folder))
  return(invisible(NULL))
}

