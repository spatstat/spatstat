#
#   pcfinhom.R
#
#   $Revision: 1.17 $   $Date: 2015/06/22 01:30:38 $
#
#   inhomogeneous pair correlation function of point pattern 
#
#

pcfinhom <- function(X, lambda=NULL, ..., r=NULL,
                     kernel="epanechnikov", bw=NULL, stoyan=0.15,
                     correction=c("translate", "Ripley"),
                     divisor=c("r","d"),
                     renormalise=TRUE,
                     normpower=1,
                     reciplambda=NULL, 
                     sigma=NULL, varcov=NULL)
{
  verifyclass(X, "ppp")
#  r.override <- !is.null(r)

  win <- X$window
  areaW <- area(win)
  npts <- npoints(X)
  
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
  
  if(is.null(bw) && kernel=="epanechnikov") {
    # Stoyan & Stoyan 1995, eq (15.16), page 285
    h <- stoyan /sqrt(npts/areaW)
    hmax <- h
    # conversion to standard deviation
    bw <- h/sqrt(5)
  } else if(is.numeric(bw)) {
    # standard deviation of kernel specified
    # upper bound on half-width
    hmax <- 3 * bw
  } else {
    # data-dependent bandwidth selection: guess upper bound on half-width
    hmax <- 2 * stoyan /sqrt(npts/areaW)
  }


  ########## intensity values #########################

  dangerous <- c("lambda", "reciplambda")
  danger <- TRUE

  if(npts == 0) {
    lambda <- reciplambda <- numeric(0)
    danger <- FALSE
  } else if(missing(lambda) && is.null(reciplambda)) {
    # No intensity data provided
    danger <- FALSE
    # Estimate density by leave-one-out kernel smoothing
    lambda <- density(X, ..., sigma=sigma, varcov=varcov,
                      at="points", leaveoneout=TRUE)
    lambda <- as.numeric(lambda)
    reciplambda <- 1/lambda
  } else if(!is.null(reciplambda)) {
    # 1/lambda values provided
    if(is.im(reciplambda)) 
      reciplambda <- safelookup(reciplambda, X)
    else if(is.function(reciplambda))
      reciplambda <- reciplambda(X$x, X$y)
    else if(is.numeric(reciplambda) && is.vector(as.numeric(reciplambda)))
      check.nvector(reciplambda, npts)
    else stop(paste(sQuote("reciplambda"),
                    "should be a vector, a pixel image, or a function"))
  } else {
    # lambda values provided
    if(is.im(lambda)) 
      lambda <- safelookup(lambda, X)
    else if(is.ppm(lambda))
      lambda <- predict(lambda, locations=X, type="trend")
    else if(is.function(lambda)) 
      lambda <- lambda(X$x, X$y)
    else if(is.numeric(lambda) && is.vector(as.numeric(lambda)))
      check.nvector(lambda, npts)
    else stop(paste(sQuote("lambda"),
                    "should be a vector, a pixel image, or a function"))
    # evaluate reciprocal
    reciplambda <- 1/lambda
  }
  
  # renormalise
  if(renormalise && npts > 0) {
    check.1.real(normpower)
    stopifnot(normpower %in% 1:2)
    renorm.factor <- (areaW/sum(reciplambda))^normpower
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

  ########## smoothing parameters for pcf ############################  
  # arguments for 'density'

  denargs <- resolve.defaults(list(kernel=kernel, bw=bw),
                              list(...),
                              list(n=length(r), from=0, to=rmax))
  
  #################################################
  
  # compute pairwise distances

  if(npts > 1) {
    close <- closepairs(X, rmax+hmax)
    dIJ <- close$d
    I <- close$i
    J <- close$j
    XI <- ppp(close$xi, close$yi, window=win, check=FALSE)
    wIJ <- reciplambda[I] * reciplambda[J]
  } else {
    undefined <- rep(NaN, length(r))
  }

  # initialise fv object
  
  df <- data.frame(r=r, theo=rep.int(1,length(r)))
  out <- fv(df, "r",
            quote(g[inhom](r)), "theo", ,
            alim,
            c("r","{%s[%s]^{pois}}(r)"),
            c("distance argument r", "theoretical Poisson %s"),
            fname=c("g", "inhom"))

  ###### compute #######

  if(any(correction=="translate")) {
    # translation correction
    if(npts > 1) {
      XJ <- ppp(close$xj, close$yj, window=win, check=FALSE)
      edgewt <- edge.Trans(XI, XJ, paired=TRUE)
      gT <- sewpcf(dIJ, edgewt * wIJ, denargs, areaW, divisor)$g
      if(renormalise) gT <- gT * renorm.factor
    } else gT <- undefined
    out <- bind.fv(out,
                   data.frame(trans=gT),
                   "{hat(%s)[%s]^{Trans}}(r)",
                   "translation-corrected estimate of %s",
                   "trans")
  }
  if(any(correction=="isotropic")) {
    # Ripley isotropic correction
    if(npts > 1) {
      edgewt <- edge.Ripley(XI, matrix(dIJ, ncol=1))
      gR <- sewpcf(dIJ, edgewt * wIJ, denargs, areaW, divisor)$g
      if(renormalise) gR <- gR * renorm.factor
    } else gR <- undefined
    out <- bind.fv(out,
                   data.frame(iso=gR),
                   "{hat(%s)[%s]^{Ripley}}(r)",
                   "isotropic-corrected estimate of %s",
                   "iso")
  }
  
  # sanity check
  if(is.null(out)) {
    warning("Nothing computed - no edge corrections chosen")
    return(NULL)
  }
  
  # which corrections have been computed?
  nama2 <- names(out)
  corrxns <- rev(nama2[nama2 != "r"])

  # default is to display them all
  formula(out) <- deparse(as.formula(paste(
                       "cbind(",
                        paste(corrxns, collapse=","),
                        ") ~ r")))
  unitname(out) <- unitname(X)
  if(danger)
    attr(out, "dangerous") <- dangerous
  return(out)
}

