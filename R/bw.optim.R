#
# bw.optim.R
#
#  Class of optimised bandwidths
#  Plotting the object displays the optimisation criterion
#
#  $Revision: 1.11 $  $Date: 2013/04/15 09:00:53 $
#

bw.optim <- function(cv, h, iopt=which.min(cv), ...,
                     cvname, hname,
                     criterion="cross-validation") {
  if(missing(cvname) || is.null(cvname)) cvname <- deparse(substitute(cv))
  if(missing(hname) || is.null(hname)) hname <- deparse(substitute(h))
  stopifnot(is.numeric(cv))
  stopifnot(is.numeric(h))
  stopifnot(length(h) == length(cv))
  result <- h[iopt]
  attr(result, "cv") <- cv
  attr(result, "h") <- h
  attr(result, "iopt") <- iopt
  attr(result, "labels") <- list(hname=hname, cvname=cvname)
  attr(result, "info") <- list(...)
  attr(result, "criterion") <- criterion
  class(result) <- "bw.optim"
  return(result)
}

print.bw.optim <- function(x, ...) {
  y <- as.numeric(x)
  names(y) <- attr(x, "labels")$hname
  print(y, ...)
  return(invisible(NULL))
}

as.data.frame.bw.optim <- function(x, ...) {
  h <- attr(x, "h")
  cv <- attr(x, "cv")
  df <- data.frame(h, cv)
  labels <- attr(x, "labels")
  colnames(df) <- labels[c("hname", "cvname")]
  return(df)
}

plot.bw.optim <- function(x, ...,
                          showopt=TRUE, optargs=list(lty=3, col="blue")) {
  xname <- short.deparse(substitute(x))
  # convert to fv object
  df <- as.data.frame(x)
  hname <- colnames(df)[1]
  cvname <- colnames(df)[2]
  yexp <- substitute(CV(h), list(CV=as.name(cvname), h=as.name(hname)))
  xfv <- fv(df,
            argu=hname,
            ylab=yexp,
            valu=cvname,
            desc=c("smoothing parameter",
                   paste(attr(x, "criterion"), "criterion")),
            fname=cvname,
            yexp=yexp)
  fvnames(xfv, ".") <- cvname
  # plot cross-validation criterion
  out <- do.call("plot.fv",
                 resolve.defaults(list(x=xfv),
                                  list(...),
                                  list(main=xname)))
  if(showopt) {
    hopt <- as.numeric(x)
    do.call("abline", append(list(v=hopt), optargs))
  }
  if(is.null(out)) return(invisible(NULL))
  return(out)
}


