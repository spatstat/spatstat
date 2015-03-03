##  vesiclesCopy.R
##  $Revision: 1.3 $ $Date: 2015/03/03 10:59:25 $

vesiclesCopy <- local({

  vfiles <- c("vesicles.txt",
              "vesicleswindow.txt",
              "vesicleswindow.csv",
              "activezone.txt")
  
  vesiclesCopy <- function(folder = getwd()) {
    oldfolder <- getwd()
    setwd(folder)
    for(fn in vfiles) {
      file.copy(
        from = system.file("rawdata", "vesicles", fn, package="spatstat"),
        to = fn,
        overwrite=TRUE)
    }
    splat("Copied files", commasep(dQuote(vfiles)), "to", dQuote(folder))
    setwd(oldfolder)
    return(invisible(NULL))
  }

  vesiclesCopy
})


