#
#   pcf.R
#
#   $Revision: 1.47 $   $Date: 2013/12/11 02:12:43 $
#
#
#   calculate pair correlation function
#   from point pattern (pcf.ppp)
#   or from estimate of K or Kcross (pcf.fv)
#   or from fasp object
#
#
pcf <- function(X, ...) {
  UseMethod("pcf")
}

pcf.ppp <- function(X, ..., r=NULL,
                    kernel="epanechnikov", bw=NULL, stoyan=0.15,
                    correction=c("translate", "Ripley"),
                    divisor=c("r", "d"))
{
  verifyclass(X, "ppp")
  r.override <- !is.null(r)

  win <- X$window
  area <- area.owin(win)
  lambda <- X$n/area
  lambda2area <- area * lambda^2

  correction.given <- !missing(correction)
  correction <- pickoption("correction", correction,
                           c(isotropic="isotropic",
                             Ripley="isotropic",
                             trans="translate",
                             translate="translate",
                             translation="translate",
                             best="best"),
                           multi=TRUE)

  correction <- implemented.for.K(correction, win$type, correction.given)

  divisor <- match.arg(divisor)
  
  # bandwidth
  if(is.null(bw) && kernel=="epanechnikov") {
    # Stoyan & Stoyan 1995, eq (15.16), page 285
    h <- stoyan /sqrt(lambda)
    hmax <- h
    # conversion to standard deviation
    bw <- h/sqrt(5)
  } else if(is.numeric(bw)) {
    # standard deviation of kernel specified
    # upper bound on half-width
    hmax <- 3 * bw
  } else {
    # data-dependent bandwidth selection: guess upper bound on half-width
    hmax <- 2 * stoyan /sqrt(lambda)
  }

  ########## r values ############################
  # handle arguments r and breaks 

  rmaxdefault <- rmax.rule("K", win, lambda)        
  breaks <- handle.r.b.args(r, NULL, win, rmaxdefault=rmaxdefault)
  if(!(breaks$even))
    stop("r values must be evenly spaced")
  # extract r values
  r <- breaks$r
  rmax <- breaks$max
  # recommended range of r values for plotting
  alim <- c(0, min(rmax, rmaxdefault))

  # arguments for 'density'
  denargs <- resolve.defaults(list(kernel=kernel, bw=bw),
                              list(...),
                              list(n=length(r), from=0, to=rmax))
  
  #################################################
  
  # compute pairwise distances
  
  close <- closepairs(X, rmax + hmax)
  dIJ <- close$d
  XI <- ppp(close$xi, close$yi, window=win, check=FALSE)

  # initialise fv object
  
  df <- data.frame(r=r, theo=rep.int(1,length(r)))
  out <- fv(df, "r",
            substitute(g(r), NULL), "theo", ,
            alim,
            c("r","%s[Pois](r)"),
            c("distance argument r", "theoretical Poisson %s"),
            fname="g")

  ###### compute #######

  if(any(correction=="translate")) {
    # translation correction
    XJ <- ppp(close$xj, close$yj, window=win, check=FALSE)
    edgewt <- edge.Trans(XI, XJ, paired=TRUE)
    gT <- sewpcf(dIJ, edgewt, denargs, lambda2area, divisor)$g
    out <- bind.fv(out,
                   data.frame(trans=gT),
                   "hat(%s)[Trans](r)",
                   "translation-corrected estimate of %s",
                   "trans")
  }
  if(any(correction=="isotropic")) {
    # Ripley isotropic correction
    edgewt <- edge.Ripley(XI, matrix(dIJ, ncol=1))
    gR <- sewpcf(dIJ, edgewt, denargs, lambda2area, divisor)$g
    out <- bind.fv(out,
                   data.frame(iso=gR),
                   "hat(%s)[Ripley](r)",
                   "isotropic-corrected estimate of %s",
                   "iso")
  }
  
  # sanity check
  if(is.null(out)) {
    warning("Nothing computed - no edge corrections chosen")
    return(NULL)
  }
  
  # default is to display all corrections
  formula(out) <- . ~ r
  #
  unitname(out) <- unitname(X)
  return(out)
}

# Smoothing Estimate of Weighted Pair Correlation
# d = vector of relevant distances
# w = vector of edge correction weights (in normal use)
# denargs = arguments to density.default
# lambda2area = constant lambda^2 * area (in normal use)

sewpcf <- function(d, w, denargs, lambda2area, divisor=c("r","d")) {
  divisor <- match.arg(divisor)
  if(divisor == "d") {
    w <- w/d
    if(!all(good <- is.finite(w))) {
      nbad <- sum(!good)
      warning(paste(nbad, "infinite or NA",
                    ngettext(nbad, "contribution was", "contributions were"),
                    "deleted from pcf estimate"))
      d <- d[good]
      w <- w[good]
    }
  }
  wtot <- sum(w)
  kden <- do.call.matched("density.default",
                  append(list(x=d, weights=w/wtot), denargs))
  r <- kden$x
  y <- kden$y * wtot
  if(divisor == "r")
    y <- y/r
  g <- y/(2 * pi * lambda2area)
  return(data.frame(r=r,g=g))
}

#
#---------- OTHER METHODS FOR pcf --------------------
#

"pcf.fasp" <- function(X, ..., method="c") {
  verifyclass(X, "fasp")
  Y <- X
  Y$title <- paste("Array of pair correlation functions",
                   if(!is.null(X$dataname)) "for",
                   X$dataname)
  # go to work on each function
  for(i in seq_along(X$fns)) {
    Xi <- X$fns[[i]]
    PCFi <- pcf.fv(Xi, ..., method=method)
    Y$fns[[i]] <- PCFi
    if(is.fv(PCFi))
      Y$default.formula[[i]] <- formula(PCFi)
  }
  return(Y)
}


pcf.fv <- local({

  callmatched <- function(fun, argue) {
    formalnames <- names(formals(fun))
    formalnames <- formalnames[formalnames != "..."]
    do.call("fun", argue[names(argue) %in% formalnames])
  }

  pcf.fv <- function(X, ..., method="c") {
    verifyclass(X, "fv")
  
    # extract r and the recommended estimate of K
    r <- with(X, .x)
    K <- with(X, .y)
    alim <- attr(X, "alim")

    # remove NA's
    ok <- !is.na(K)
    K <- K[ok]
    r <- r[ok]
    switch(method,
           a = {
             ss <- callmatched(smooth.spline,
                               list(x=r, y=K, ...))
             dK <- predict(ss, r, deriv=1)$y
             g <- dK/(2 * pi * r)
           },
           b = {
             y <- K/(2 * pi * r)
             y[!is.finite(y)] <- 0
             ss <- callmatched(smooth.spline,
                               list(x=r, y=y, ...))
             dy <- predict(ss, r, deriv=1)$y
             g <- dy + y/r
           },
           c = {
             z <- K/(pi * r^2)
             z[!is.finite(z)] <- 1
             ss <- callmatched(smooth.spline,
                               list(x=r, y=z, ...))
             dz <- predict(ss, r, deriv=1)$y
             g <- (r/2) * dz + z
           },
           d = {
             z <- sqrt(K)
             z[!is.finite(z)] <- 0
             ss <- callmatched(smooth.spline,
                               list(x=r, y=z, ...))
             dz <- predict(ss, r, deriv=1)$y
             g <- z * dz/(pi * r)
           },
           stop(paste("unrecognised method", sQuote(method)))
           )

    # pack result into "fv" data frame
    Z <- fv(data.frame(r=r,
                       theo=rep.int(1, length(r)),
                       pcf=g),
            "r", substitute(g(r), NULL), "pcf", . ~ r, alim,
            c("r", "%s[pois](r)", "%s(r)"),
            c("distance argument r",
              "theoretical Poisson value of %s",
              "estimate of %s by numerical differentiation"),
            fname="g")
    unitname(Z) <- unitname(X)
    return(Z)
  }

  pcf.fv
})

