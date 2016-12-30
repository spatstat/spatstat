#'
#'     timed.R
#'
#'   Timed objects
#'
#'   $Revision: 1.1 $ $Date: 2016/12/30 03:24:47 $

timed <- function(x, ..., starttime=NULL, timetaken=NULL) {
  if(is.null(starttime) && is.null(timetaken)) # time starts now.
    starttime <- proc.time()
  # evaluate expression if any
  object <- x
  if(is.null(timetaken))
    timetaken <- proc.time() - starttime
  if(!inherits(object, "timed"))
    class(object) <- c("timed", class(object))
  attr(object, "timetaken") <- timetaken
  return(object)
}

print.timed <- function(x, ...) {
  # strip the timing information and print the rest.
  taken <- attr(x, "timetaken")
  cx <- class(x)
  attr(x, "timetaken") <- NULL
  class(x) <- cx[cx != "timed"]
  NextMethod("print")
  # Now print the timing info
  cat(paste("\nTime taken:", codetime(taken), "\n"))
  return(invisible(NULL))
}

timeTaken <- function(..., warn=TRUE) {
  allargs <- list(...)
  hastime <- sapply(allargs, inherits, what="timed")
  if(warn && !all(hastime))
    warning("Some arguments did not contain timing information")
  times <- sapply(allargs[hastime], attr, which="timetaken")
  tottime <- rowSums(times)
  class(tottime) <- "proc_time"
  return(tottime)
}

#'  ..............  codetime ....................................
#'  Basic utility for converting times in seconds to text strings

codetime <- local({
  uname <- c("min", "hours", "days", "years",
             "thousand years", "million years", "billion years")
  u1name <- c("min", "hour", "day", "year",
             "thousand years", "million years", "billion years")
  multiple <- c(60, 60, 24, 365, 1e3, 1e3, 1e3)
  codehms <- function(x) {
    sgn <- if(x < 0) "-" else ""
    x <- round(abs(x))
    hours <- x %/% 3600
    mins  <- (x %/% 60) %% 60
    secs  <- x %% 60
    h <- if(hours > 0) paste(hours, ":", sep="") else ""
    started <- (hours > 0)
    m <- if(mins > 0) {
      paste(if(mins < 10 && started) "0" else "", mins, ":", sep="")
    } else if(started) "00:" else ""
    started <- started | (mins > 0)
    s <- if(secs > 0) {
      paste(if(secs < 10 && started) "0" else "", secs, sep="")
    } else if(started) "00" else "0"
    if(!started) s <- paste(s, "sec")
    paste(sgn, h, m, s, sep="")
  }
  codetime <- function(x, hms=TRUE, what=c("elapsed","user","system")) {
    if(inherits(x, "proc_time")) x <- summary(x)[[match.arg(what)]] 
    if(!is.numeric(x) || length(x) != 1)
      stop("codetime: x must be a proc_time object or a single number")
    sgn <- if(x < 0) "-" else ""
    x <- abs(x)
    if(x < 60)
      return(paste(sgn, signif(x, 3), " sec", sep=""))
    # more than 1 minute: round to whole number of seconds
    x <- round(x)
    if(hms && (x < 60 * 60 * 24))
      return(paste(sgn, codehms(x), sep=""))
    u <- u1 <- "sec"
    for(k in seq_along(multiple)) {
      if(x >= multiple[k]) {
        x <- x/multiple[k]
        u <- uname[k]
        u1 <- u1name[k]
      } else break
    }
    xx <- round(x, 1)
    ux <- if(xx == 1) u1 else u
    paste(sgn, xx, " ", ux, sep="")
  }
  codetime
})

