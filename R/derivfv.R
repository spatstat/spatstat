#
# derivfv.R
#
# differentiation for fv objects
#
#  $Revision: 1.6 $ $Date: 2014/10/24 00:22:30 $
#

deriv.fv <- local({

  derivative <- function(y, r, ...) {
    ss <- smooth.spline(r, y, ...)
    predict(ss, r, deriv=1)$y
  }
  
  deriv.fv <- function(expr, which="*", ...,
                       method=c("spline", "numeric"),
                       kinks=NULL,
                       periodic=FALSE,
                       Dperiodic=periodic) {
    f <- expr
    method <- match.arg(method)

    ## select columns
    ##  if(length(which) == 1 && which %in% .Spatstat.FvAbbrev) {
    if(length(which) == 1) {
      if(which == ".x")
        stop("Cannot smooth the function argument")
      which <- fvnames(f, which)
    }
    
    if(any(nbg <- !(which %in% names(f)))) 
      stop(paste("Unrecognised column",
                 ngettext(sum(nbg), "name", "names"),
                 commasep(sQuote(which[nbg])), 
                 "in argument", sQuote("which")))
    relevant <- names(f) %in% which
    ## get 
    rname <- fvnames(f, ".x")
    df <- as.data.frame(f)
    rpos <- which(colnames(df) == rname)
    rvals <- df[,rpos]
    yvals <- df[,relevant,drop=FALSE]
    nr <- length(rvals)
    ##
    if(Dperiodic) {
      ## Derivative should be periodic
      ## Recycle data to imitate periodicity
      DR <- diff(range(rvals))
      rvals <- c(rvals[-nr] - DR, rvals, rvals[-1] + DR)
      yleft <- yvals[-nr, , drop=FALSE]
      yright <-  yvals[-1, , drop=FALSE]
      if(!periodic) {
        ## original data are not periodic (e.g. cdf of angular variable)
        ## but derivative must be periodic
        jump <- matrix(as.numeric(yvals[nr,] - yvals[1, ]),
                       nr-1, ncol(yvals), byrow=TRUE)
        yleft <- yleft - jump
        yright <- yright + jump
      }
      yvals <- rbind(yleft, yvals, yright)
      actual <- nr:(2*nr - 1)
      NR <- length(rvals)
    } else {
      NR <- nr
      actual <- 1:nr
    }
    ## cut x axis into intervals?
    if(is.null(kinks)) {
      cutx <- factor(rep(1, NR))
    } else {
      rr <- range(rvals)
      if(periodic) 
        kinks <- c(kinks-DR, kinks, kinks+DR)
      breaks <- sort(unique(kinks))
      if(breaks[1] > rr[1]) breaks <- c(rr[1], breaks)
      if(max(breaks) < rr[2]) breaks <- c(breaks, rr[2])
      cutx <- cut(rvals, breaks=breaks, include.lowest=TRUE)
    }
    ## process
    for(segment in levels(cutx)) {
      ii <- (cutx == segment)
      yy <- yvals[ii, , drop=FALSE]
      switch(method,
             numeric = {
               dydx <- apply(yy, 2, diff)/diff(rvals[ii])
               nd <- nrow(dydx)
               dydx <- rbind(dydx, dydx[nd, ])
             },
             spline = {
               dydx <- apply(yy, 2, derivative, 
                             r=rvals[ii], ...)
         })
      df[ii[actual], relevant] <- dydx[ actual, ]
    }
    ## pack up
    result <- f
    result[,] <- df
    ## tweak name of function
    if(!is.null(yl <- attr(f, "ylab")))
      attr(result, "ylab") <- substitute(bold(D)~Fx, list(Fx=yl))
    if(!is.null(ye <- attr(f, "yexp")))
      attr(result, "yexp") <- substitute(bold(D)~Fx, list(Fx=ye))
    ## tweak mathematical labels
    attr(result, "labl")[relevant]  <-
      paste0("bold(D)~", attr(f, "labl")[relevant])
    return(result)
  }

  deriv.fv
})


increment.fv <- function(f, delta) {
  stopifnot(is.fv(f))
  check.1.real(delta)
  stopifnot(delta > 0)
  half <- delta/2
  xx <- with(f, .x)
  ynames <- fvnames(f, ".")
  yy <- as.data.frame(lapply(ynames,
                             function(a, xx, f, h) {
                               g <- as.function(f, value=a)
                               g(xx+h)-g(xx-h)
                             },
                             xx=xx, f=f, h=half))
  Y <- f
  Y[,ynames] <- yy
  ## tweak name of function
  if(!is.null(yl <- attr(f, "ylab")))
    attr(Y, "ylab") <- substitute(Delta~Fx, list(Fx=yl))
  if(!is.null(ye <- attr(f, "yexp")))
    attr(Y, "yexp") <- substitute(Delta~Fx, list(Fx=ye))
  ## tweak mathematical labels
  relevant <- colnames(Y) %in% ynames
  attr(Y, "labl")[relevant]  <-
      paste0("Delta~", attr(f, "labl")[relevant])
  ## tweak recommended range
  attr(Y, "alim") <- intersect.ranges(attr(f, "alim"),
                                      range(xx) + c(1,-1)*half)
  return(Y)
}

  
