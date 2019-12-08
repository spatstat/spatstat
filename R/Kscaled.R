#
#	Kscaled.R	Estimation of K function for locally-scaled process
#
#	$Revision: 1.17 $	$Date: 2019/12/08 04:29:28 $
#

"Lscaled" <- function(...) {
  K <- Kscaled(...)
  L <- eval.fv(sqrt(pmax.int(K,0)/pi))
  # relabel the fv object
  L <- rebadge.fv(L, quote(L[scaled](r)), c("L","scaled"))
  attr(L, "labl") <- attr(K, "labl")
  return(L)  
}

"Kscaled"<-
  function (X, lambda=NULL, ..., r = NULL, breaks = NULL,
            rmax = 2.5,
            correction=c("border", "isotropic", "translate"),
            renormalise=FALSE, normpower=1,
            sigma=NULL, varcov=NULL)
{
  verifyclass(X, "ppp")
#  rfixed <- !missing(r) || !missing(breaks)

  ## determine basic parameters
  W <- X$window
  npts <- X$n
  areaW <- area(W)
  halfdiameter <- diameter(W)/2
  
  ## match corrections
  correction.given <- !missing(correction) && !is.null(correction)
  correction <- pickoption("correction", correction,
                           c(none="none",
                             border="border",
                             isotropic="isotropic",
                             Ripley="isotropic",
                             trans="translate",
                             translate="translate",
                             translation="translate",
                             best="best"),
                           multi=TRUE)

#  best.wanted <- ("best" %in% correction)
  correction <- implemented.for.K(correction, W$type, correction.given)

  ###########################################################
  ## DETERMINE WEIGHTS AND VALIDATE
  ##

  if(missing(lambda)) {
    ## No intensity data provided
    ## Estimate density by leave-one-out kernel smoothing
    lambda <- density(X, ..., sigma=sigma, varcov=varcov,
                      at="points", leaveoneout=TRUE)
    lambda <- as.numeric(lambda)
  } else {
    ## lambda values provided
    if(is.im(lambda)) 
      lambda <- safelookup(lambda, X)
    else if(is.function(lambda)) 
      lambda <- lambda(X$x, X$y)
    else if(is.ppm(lambda)) 
      lambda <- safelookup(predict(lambda, type="trend"), X)
    else if(!is.numeric(lambda) || !is.null(dim(lambda)))
      stop(paste(sQuote("lambda"),
                 "should be a vector, a pixel image, a function or a ppm"))
    check.nvector(lambda, npts)
  }

  if(renormalise) {
    ## renormalise. Here we only need half the power ;-)
    check.1.real(normpower)
    stopifnot(normpower %in% 1:2) 
    renorm.factor <- (areaW/sum(1/lambda))^(normpower/2)
    lambda <- lambda/renorm.factor
  }     
  ## Calculate range of r values using max lambda
  sra <- sqrt(range(lambda))
  minrescale <- sra[1]
  maxrescale <- sra[2]

  ## convert arguments to absolute distances 
  absr <- if(!is.null(r)) r/maxrescale else NULL
  absrmaxdefault <- min(rmax.rule("K", W), rmax/maxrescale)
  absbreaks <-
    if(!is.null(breaks)) scalardilate(breaks, 1/maxrescale) else NULL
  ## determine absolute distances
  absbreaks <- handle.r.b.args(absr, absbreaks, W, rmaxdefault=absrmaxdefault)
  absr <- absbreaks$r
  ## convert to rescaled distances
  breaks <- scalardilate(absbreaks, maxrescale)
  r <- breaks$r
  rmax <- breaks$max
  ## recommended range of scaled r values
  alim <- c(0, min(rmax, maxrescale * absrmaxdefault))
  rthresh <- minrescale * halfdiameter
  ## maximum absolute distance ever needed
  maxabsdist <- min(rmax/minrescale, halfdiameter)
  
  ## this will be the output data frame
  K <- data.frame(r=r, theo= pi * r^2)
  desc <- c("distance argument r", "theoretical Poisson %s")
  K <- fv(K, "r", quote(K[scaled](r)),
          "theo", , alim,
          c("r","{%s[%s]^{pois}}(r)"),
          desc,
          fname=c("K", "scaled"))
        
  ## identify all relevant close pairs
  needXI <- any(correction %in% c("translate", "isotropic"))
  close <- closepairs(X, maxabsdist, what=if(needXI) "all" else "ijd")
  I <- close$i
  J <- close$j
  ## locally-scaled distances
  sqrtLambda <- sqrt(lambda)
  lamIJ <- (sqrtLambda[I] + sqrtLambda[J])/2
  absDIJ <- close$d
  DIJ <- absDIJ * lamIJ
  ## first point of each pair
  XI <- if(needXI) ppp(close$xi, close$yi, window=W, check=FALSE) else NULL
  
  if(any(correction == "none")) {
    ## uncorrected! For demonstration purposes only!
    wh <- whist(DIJ, breaks$val)  # no weights
    Kun <- cumsum(wh)/npts
    K <- bind.fv(K, data.frame(un=Kun), "{hat(%s)[%s]^{un}}(r)",
                 "uncorrected estimate of %s",
                 "un")
  }
  
  if(any(correction == "border")) {
    ## border method
    ## Compute SCALED distances to boundary
    b <- bdist.points(X) * sqrtLambda
    bI <- b[I]
    ## apply reduced sample algorithm to scaled distances
    RS <- Kount(DIJ, bI, b, breaks)
    Kb <- RS$numerator/RS$denom.count
    Kb[r > rthresh] <- NA
    K <- bind.fv(K, data.frame(border=Kb), "{hat(%s)[%s]^{bord}}(r)",
                 "border-corrected estimate of %s",
                 "border")
  }

  if(any(correction == "translate")) {
    ## translation correction
    XJ <- ppp(close$xj, close$yj, window=W, check=FALSE)
    edgewt <- edge.Trans(XI, XJ, paired=TRUE)
    wh <- whist(DIJ, breaks$val, edgewt)
    Ktrans <- cumsum(wh)/npts
    Ktrans[r >= rthresh] <- NA
    K <- bind.fv(K, data.frame(trans=Ktrans), "{hat(%s)[%s]^{trans}}(r)",
                 "translation-corrected estimate of %s",
                 "trans")
  }
  if(any(correction == "isotropic")) {
    ## Ripley isotropic correction (using UN-SCALED distances)
    edgewt <- edge.Ripley(XI, matrix(absDIJ, ncol=1))
    wh <- whist(DIJ, breaks$val, edgewt)
    Kiso <- cumsum(wh)/npts
    Kiso[r >= rthresh] <- NA
    K <- bind.fv(K, data.frame(iso=Kiso), "{hat(%s)[%s]^{iso}}(r)",
                 "Ripley isotropic correction estimate of %s",
                 "iso")
  }
  ## default plot will display all edge corrections
  formula(K) <- . ~ r
  nama <- rev(colnames(K))
  fvnames(K, ".") <- nama[!(nama %in% c("r", "rip", "ls"))]
  ##
  unitname(K) <- c("normalised unit", "normalised units")
  return(K)
}
	
