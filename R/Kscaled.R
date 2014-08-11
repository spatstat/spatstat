#
#	Kscaled.R	Estimation of K function for locally-scaled process
#
#	$Revision: 1.7 $	$Date: 2013/02/07 09:58:14 $
#

"Lscaled" <- function(...) {
  K <- Kscaled(...)
  L <- eval.fv(sqrt(pmax(K,0)/pi))
  # relabel the fv object
  L <- rebadge.fv(L, substitute(Lscaled(r), NULL), "Lscaled")
  return(L)  
}

"Kscaled"<-
  function (X, lambda=NULL, ..., r = NULL, breaks = NULL, 
         correction=c("border", "isotropic", "translate"),
            sigma=NULL, varcov=NULL)
{
    verifyclass(X, "ppp")
    rfixed <- !missing(r) || !missing(breaks)

    # determine basic parameters
    W <- X$window
    npts <- X$n
    area <- area.owin(W)

    rmaxdefault <- rmax.rule("K", W, npts/area) * sqrt(npts/area)
    breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
    r <- breaks$r
    rmax <- breaks$max

    # match corrections
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

    best.wanted <- ("best" %in% correction)
    correction <- implemented.for.K(correction, W$type, correction.given)

    ###########################################################
    # DETERMINE WEIGHTS AND VALIDATE
    #

    if(missing(lambda)) {
      # No intensity data provided
      # Estimate density by leave-one-out kernel smoothing
      lambda <- density(X, ..., sigma=sigma, varcov=varcov,
                        at="points", leaveoneout=TRUE)
      lambda <- as.numeric(lambda)
    } else {
        # lambda values provided
      if(is.im(lambda)) 
        lambda <- safelookup(lambda, X)
      else if(is.function(lambda)) 
        lambda <- lambda(X$x, X$y)
      else if(is.numeric(lambda) && is.vector(as.numeric(lambda)))
        check.nvector(lambda, npts)
      else stop(paste(sQuote("lambda"),
                      "should be a vector, a pixel image, or a function"))
    }
    
    # recommended range of r values
    alim <- c(0, min(rmax, rmaxdefault))
        
    # this will be the output data frame
    K <- data.frame(r=r, theo= pi * r^2)
    desc <- c("distance argument r", "theoretical Poisson %s")
    K <- fv(K, "r", substitute(Kscaled(r), NULL),
            "theo", , alim, c("r","%s[pois](r)"), desc, fname="Kscaled")
        
    # identify all close pairs
    rmax <- max(r)
    close <- closepairs(X, rmax)
    I <- close$i
    J <- close$j
    # locally-scaled distances
    lamIJ <- (sqrt(lambda[I]) + sqrt(lambda[J]))/2
    absDIJ <- close$d
    DIJ <- absDIJ * lamIJ

    XI <- ppp(close$xi, close$yi, window=W, check=FALSE)
  
  if(any(correction == "none")) {
    # uncorrected! For demonstration purposes only!
    wh <- whist(DIJ, breaks$val)  # no weights
    Kun <- cumsum(wh)/npts
    K <- bind.fv(K, data.frame(un=Kun), "%s[un](r)",
                 "uncorrected estimate of %s",
                 "un")
  }
  
  if(any(correction == "border")) {
  # border method
  # Compute SCALED distances to boundary
    b <- bdist.points(X) * sqrt(lambda)
    I <- close$i
    bI <- b[I]
  # apply reduced sample algorithm
    RS <- Kount(DIJ, bI, b, breaks)
    if(any(correction == "border")) {
      Kb <- RS$numerator/RS$denom.count
      K <- bind.fv(K, data.frame(border=Kb), "%s[bord](r)",
                   "border-corrected estimate of %s",
                   "border")
    }
  }

  if(any(correction == "translate")) {
    # translation correction
    XJ <- ppp(close$xj, close$yj, window=W, check=FALSE)
    edgewt <- edge.Trans(XI, XJ, paired=TRUE)
    wh <- whist(DIJ, breaks$val, edgewt)
    Ktrans <- cumsum(wh)/npts
    h <- diameter(W)/2
    Ktrans[r >= h] <- NA
    K <- bind.fv(K, data.frame(trans=Ktrans), "%s[trans](r)",
                 "translation-corrected estimate of %s",
                 "trans")
  }
  if(any(correction == "isotropic")) {
    # Ripley isotropic correction (using UN-SCALED distances)
    edgewt <- edge.Ripley(XI, matrix(absDIJ, ncol=1))
    wh <- whist(DIJ, breaks$val, edgewt)
    Kiso <- cumsum(wh)/npts
    h <- diameter(W)/2
    Kiso[r >= h] <- NA
    K <- bind.fv(K, data.frame(iso=Kiso), "%s[iso](r)",
                 "Ripley isotropic correction estimate of %s",
                 "iso")
  }
  # default plot will display all edge corrections
  formula(K) <- . ~ r
  nama <- rev(colnames(K))
  fvnames(K, ".") <- nama[!(nama %in% c("r", "rip", "ls"))]
  #
  unitname(K) <- c("normalised unit", "normalised units")
  return(K)
}
	
