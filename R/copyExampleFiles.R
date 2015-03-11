##  copyExampleFiles.R
##  $Revision: 1.10 $ $Date: 2015/03/11 05:58:50 $

copyExampleFiles <- function(which, folder=getwd()) {
  choices <- dir(system.file("rawdata", package="spatstat"))
  if(missing(which) || is.null(which)) {
    message(paste("Choices are: which=", commasep(sQuote(choices), " or ")))
    return(invisible(NULL))
  }
  if(!interactive())
    stop("Copying files requires an interactive session (by CRAN Policies).")
  whichdata <- match.arg(which, choices)
  sourcefolder <- system.file("rawdata", whichdata, package="spatstat")
  sourcefiles <- dir(sourcefolder)
  if(length(sourcefiles) == 0)
      stop("No files available")
  # set directory
  oldfolder <- getwd()
  setwd(folder)
  on.exit(setwd(oldfolder))
  # Warn user:
  foldername <- if(identical(folder, oldfolder)) "the current folder" else
                 paste("the folder", dQuote(folder))
  splat("You are about to have been copying", 
        ngettext(length(sourcefiles), "file", "files"),
        commasep(dQuote(sourcefiles)), "to",
        paste0(foldername, "."),
        "This may overwrite existing files.")
  # Ask permission:
  answer <- readline("Do you want to continue? (y/n)[y] ")
  if(!tolower(substr(answer, 1, 1)) %in% c("", "y")) {
    splat("Aborting...")
    return(invisible(NULL))
  }
  # 
  for(fn in sourcefiles) {
    frompath <- file.path(sourcefolder, fn)
    file.copy(from = frompath, to = fn, overwrite=TRUE)
  }
  splat("Copying completed.")
  return(invisible(NULL))
}
