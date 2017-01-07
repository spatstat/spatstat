#'
#'        bugtable.R
#' 
#'    $Revision: 1.3 $ $Date: 2017/01/07 04:20:31 $

bugfixes <- function(sinceversion=NULL, sincedate=NULL,
                     package="spatstat",
                     show=TRUE) {
  if(!is.null(sincedate) && package != "spatstat") {
    #' news items after specified date
    ne <- news(package=package)
    if(is.null(ne) || is.null(ne$Date) || anyNA(ne$Date))
      stop(paste(if(is.null(ne)) "News" else "Date",
                    "information is not available for package",
                 sQuote(package)), call.=FALSE)
    a <- eval(substitute(news(Date >= SD & grepl("^BUG", Category),
                              package=package),
                         list(SD=sincedate)))
  } else {
    #' determine a corresponding version number
    if(is.null(sinceversion) && is.null(sincedate)) {
      #' default is latest version
      dfile <- system.file("DESCRIPTION", package=package)
      sinceversion <- read.dcf(file=dfile, fields="Version")
    } else if(!is.null(sincedate) && package == "spatstat") {
      #' read spatstat release history table
      fname <- system.file("doc", "packagesizes.txt", package="spatstat")
      p <- read.table(fname, header=TRUE, stringsAsFactors=FALSE)
      #' find earliest package version on or after the given date
      imin <- with(p, min(which(as.Date(date) >= sincedate)))
      sinceversion <- p[imin, "version"]
    }
    a <- eval(substitute(news(Version >= sv & grepl("^BUG", Category),
                              package=package),
                         list(sv=sinceversion)))
  }
  if(!is.data.frame(a)) return(NULL)
  #' split each entry into lines
  alines <- strsplit(a$Text, "\n")
  #' extract first line
  f <- unname(sapply(alines, "[", i=1L))
  #' extract body
  b <- unname(lapply(alines, "[", i=-1L))
  b <- unname(sapply(b, paste, collapse="\n"))
  #' extract header from first line
  h <- unname(sapply(strsplit(f, ":"), "[", i=1L))
  h <- unname(sapply(strsplit(h, ","), "[", i=1L))
  h <- unname(sapply(strsplit(h, " "), "[", i=1L))
  #' sort by header
  oo <- order(h, f)
  #' rebuild
  z <- data.frame(Header=h[oo],
                  Firstline=f[oo],
                  Body=b[oo],
                  Version=a$Version[oo],
                  stringsAsFactors=FALSE)
  class(z) <- c("bugtable", class(z))
  if(!show) return(z)
  page(z, method="print")
  return(invisible(z))
}

class(bugfixes) <- "autoexec"
              
print.bugtable <- function(x, ...) {
  hprev <- ""
  for(i in seq_len(nrow(x))) {
    h <- x$Header[i]
    f <- x$Firstline[i]
    if(h != hprev) {
      # new main header
      cat("\n***", h, "***\n", fill=TRUE)
    }
    cat(x$Version[i], ":", f, fill=TRUE)
    cat(x$Body[i], "\n", fill=TRUE)
    hprev <- h
  }
  return(invisible(NULL))
}
