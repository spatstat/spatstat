#'
#'        bugtable.R
#' 
#'    $Revision: 1.16 $ $Date: 2025/03/21 07:17:48 $

bugfixes <- function(sinceversion=NULL, sincedate=NULL,
                     package=spatstat.family(),
                     show=TRUE) {
  package <- as.character(package)
  is.spat <- package %in% spatstat.family(TRUE, TRUE)
  if("book" %in% list(sinceversion, sincedate)) {
    ## special usage
    sinceversion <- NULL
    sincedate    <- "2015-06-05"
  }
  n <- length(package)
  z <- NULL
  for(i in seq_len(n)) {
    #' extract bug reports for package [i]
    packi <- package[i]
    if("all" %in% list(sinceversion, sincedate)) {
      ## no date constraint - show all news
      a <- eval(substitute(news(grepl("^BUG", Category),
                                package=packi)))
    } else if(!is.null(sincedate) && !is.spat[i]) {
      #' news items after specified date
      #' non-spatstat package
      ne <- news(package=packi)
      if(is.null(ne) || is.null(ne$Date) || anyNA(ne$Date))
        stop(paste(if(is.null(ne)) "News" else "Date",
                   "information is not available for package",
                   sQuote(packi)), call.=FALSE)
      a <- eval(substitute(news(Date >= SD & grepl("^BUG", Category),
                                package=packi),
                           list(SD=sincedate)))
    } else {
      #' Determine a starting version number (for this package i)
      #' using the arguments 'sincedate' and 'sinceversion'
      #' (whichever is more restrictive)
      startversion <- NULL
      if(is.null(sinceversion) && is.null(sincedate)) {
        #' default is latest version
        dfile <- system.file("DESCRIPTION", package=packi)
        startversion <- read.dcf(file=dfile, fields="Version")
      } else if(!is.null(sincedate) && is.spat[i]) {
        #' news items on or after specified date 'sincedate'
        #' read release history table
        fname <- system.file("doc", "packagesizes.txt", package=packi)
        if(nchar(fname) > 0) {
          ## File was found.
          p <- try(read.table(fname, header=TRUE, stringsAsFactors=FALSE,
                              col.names=c("date", "version",
                                          "nhelpfiles", "nobjects", "ndatasets",
                                          "Rlines", "srclines")))
          if(!inherits(p, "try-error") && nrow(p) > 0) {
            #' Table data were read in successfully and are non-empty.
            #' Find entries on or after the given date
            rowsAfter <- with(p, which(as.Date(date) >= sincedate))
            if(length(rowsAfter) > 0) {
              #' Find earliest entry 
              imin <- min(rowsAfter)
              #' Start from this version
              startversion <- p[imin, "version"]
            }
          }
        }
      }
      if(!is.null(sinceversion)) {
        ## Argument 'sinceversion' was given
        if(is.null(startversion) ||
           as.package_version(sinceversion) > startversion) {
          ## 'sinceversion' overrides calculated 'startversion'
          ## for this package
          startversion <- sinceversion
        }
      }
      if(is.null(startversion)) {
        #' There is no relevant package history
        a <- NULL
      } else {
        a <- eval(substitute(news(Version >= sv & grepl("^BUG", Category),
                                  package=packi),
                             list(sv=startversion)))
      }
    }
    #' convert format
    if(is.data.frame(a) && nrow(a) > 0) {
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
      #' rebuild
      zi <- data.frame(Header=h,
                       Firstline=f,
                       Body=b,
                       Version=a$Version,
                       stringsAsFactors=FALSE)
      if(n > 1) zi$Package <- packi
      z <- if(is.null(z)) zi else rbind(z, zi)
    } 
  }
  if(is.null(z)) return(NULL)
  #' sort by header
  oo <- with(z, order(Header, Firstline))
  z <- z[oo, ]
  #' wrap up
  class(z) <- c("bugtable", class(z))
  if(show) page(z, method="print")
  return(invisible(z))
}

class(bugfixes) <- "autoexec"
              
print.bugtable <- function(x, ...) {
  hprev <- ""
  haspack <- "Package" %in% colnames(x)
  for(i in seq_len(nrow(x))) {
    h <- x$Header[i]
    f <- x$Firstline[i]
    v <- x$Version[i]
    p <- if(haspack) x$Package[i] else NULL
    b <- x$Body[i]
    if(h != hprev) {
      # new main header
      cat("\n***", h, "***\n", fill=TRUE)
    }
    if(haspack) cat(p, v, ":", f, fill=TRUE) else cat(v, ":", f, fill=TRUE)
    cat(b, "\n", fill=TRUE)
    hprev <- h
  }
  return(invisible(NULL))
}
