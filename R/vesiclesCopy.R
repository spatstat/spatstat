##  vesiclesCopy.R
##  $Revision: 1.4 $ $Date: 2015/03/03 11:07:04 $

vesiclesCopy <- local({

  vfiles <- c("vesicles.txt",
              "vesicleswindow.txt",
              "vesicleswindow.csv",
              "activezone.txt")
  
  vesiclesCopy <- function(folder = getwd()) {
    oldfolder <- getwd()
    setwd(folder)
    on.exit(setwd(oldfolder))
    for(fn in vfiles) {
      file.copy(
        from = system.file("rawdata", "vesicles", fn, package="spatstat"),
        to = fn,
        overwrite=TRUE)
    }
    splat("Copied files", commasep(dQuote(vfiles)), "to", dQuote(folder))
    return(invisible(NULL))
  }

  vesiclesCopy
})


