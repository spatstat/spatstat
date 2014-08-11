#
#	Kinhom.S	Estimation of K function for inhomogeneous patterns
#
#	$Revision: 1.63 $	$Date: 2013/02/07 09:58:14 $
#
#	Kinhom()	compute estimate of K_inhom
#
#
#       Reference:
#            Non- and semiparametric estimation of interaction
#	     in inhomogeneous point patterns
#            A.Baddeley, J.Moller, R.Waagepetersen
#            Statistica Neerlandica 54 (2000) 329--350.
#
# -------- functions ----------------------------------------
#	Kinhom()	compute estimate of K
#                       using various edge corrections
#
#       Kwtsum()         internal routine for border correction
#
# -------- standard arguments ------------------------------	
#	X		point pattern (of class 'ppp')
#
#	r		distance values at which to compute K	
#
#       lambda          vector of intensity values for points of X
#
# -------- standard output ------------------------------
#      A data frame (class "fv") with columns named
#
#	r:		same as input
#
#	trans:		K function estimated by translation correction
#
#	iso:		K function estimated by Ripley isotropic correction
#
#	theo:		K function for Poisson ( = pi * r ^2 )
#
#	border:		K function estimated by border method
#			(denominator = sum of weights of points)
#
#       bord.modif:	K function estimated by border method
#			(denominator = area of eroded window)
#
# ------------------------------------------------------------------------

"Linhom" <- function(...) {
  K <- Kinhom(...)
  L <- eval.fv(sqrt(pmax(K,0)/pi))
  # relabel the fv object
  L <- rebadge.fv(L, quote(Linhom(r)), "Linhom",
                  names(K), new.labl=attr(K, "labl"))
  #
  return(L)  
}

"Kinhom"<-
  function (X, lambda=NULL, ..., r = NULL, breaks = NULL, 
         correction=c("border", "bord.modif", "isotropic", "translate"),
            renormalise=TRUE,
            normpower=1,
            nlarge = 1000, 
            lambda2=NULL,
            reciplambda=NULL, reciplambda2=NULL,
            sigma=NULL, varcov=NULL)
{
    verifyclass(X, "ppp")
    nlarge.given <- !missing(nlarge)
    rfixed <- !missing(r) || !missing(breaks)

    # determine basic parameters
    W <- X$window
    npts <- npoints(X)
    area <- area.owin(W)

    rmaxdefault <- rmax.rule("K", W, npts/area)
    breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
    r <- breaks$r
    rmax <- breaks$max

    # match corrections
    correction.given <- !missing(correction) && !is.null(correction)
    correction <- pickoption("correction", correction,
                             c(none="none",
                               border="border",
                               "bord.modif"="bord.modif",
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
    # The matrix 'lambda2' or 'reciplambda2' is sufficient information
    # unless we want the border correction.
    lambda2.given    <- !is.null(lambda2) || !is.null(reciplambda2)
    lambda2.suffices <- !any(correction %in% c("bord", "bord.modif"))

    # Use matrix of weights if it was provided and if it is sufficient
    if(lambda2.suffices && lambda2.given) {
      if(!is.null(reciplambda2)) 
        check.nmatrix(reciplambda2, npts)
      else {
        check.nmatrix(lambda2, npts)
        reciplambda2 <- 1/lambda2
      }
      # renormalise
      if(renormalise) {
        check.1.real(normpower)
        stopifnot(normpower %in% 1:2)
        renorm.factor <- (area^2/sum(reciplambda2))^(normpower/2)
      } 
    } else {
      # Vector lambda or reciplambda is required
      if(missing(lambda) && is.null(reciplambda)) {
        # No intensity data provided
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
      if(renormalise) {
        check.1.real(normpower)
        stopifnot(normpower %in% 1:2)
        renorm.factor <- (area/sum(reciplambda))^normpower
      } 
    }

    
    # recommended range of r values
    alim <- c(0, min(rmax, rmaxdefault))
        
  ###########################################
  # Efficient code for border method
  # Usable only if r values are evenly spaced from 0 to rmax
  # Invoked automatically if number of points is large

    can.do.fast <- breaks$even  && missing(lambda2)
    borderonly <- all(correction == "border" | correction == "bord.modif")
    large.n    <- (npts >= nlarge)
    demand.best <- correction.given && best.wanted
    large.n.trigger <- large.n && !demand.best
    will.do.fast <- can.do.fast && (borderonly || large.n.trigger)
    asked      <- borderonly || (nlarge.given && large.n.trigger)

    if(will.do.fast && !asked)
      message(paste("number of data points exceeds",
                    nlarge, "- computing border estimate only"))

    if(asked && !can.do.fast) {
      whynot <-
        if(!(breaks$even)) "r values not evenly spaced" else
        if(!missing(lambda)) "matrix lambda2 was given" else NULL
      warning(paste(c("cannot use efficient code", whynot), sep="; "))
    }
    
    if(will.do.fast) {
    # restrict r values to recommended range, unless specifically requested
      if(!rfixed) 
        r <- seq(from=0, to=alim[2], length.out=length(r))
      K <- Kborder.engine(X, max(r), length(r), correction, reciplambda)
      # tweak labels
      K <- rebadge.fv(K, substitute(K[inhom](r), NULL), "K[inhom]")
      K <- tweak.fv.entry(K, "theo", new.labl="{%s^{pois}}(r)")
      K <- tweak.fv.entry(K, "border", new.labl="hat(%s^{bord})(r)")
      K <- tweak.fv.entry(K, "bord.modif", new.labl="hat(%s^{bordm})(r)")
      return(K)
    }

  ###########################################
  # Slower code
  ###########################################
        
        
    # this will be the output data frame
    K <- data.frame(r=r, theo= pi * r^2)
    desc <- c("distance argument r", "theoretical Poisson %s")
    K <- fv(K, "r", substitute(K[inhom](r), NULL),
            "theo", , alim, c("r","{%s^{pois}}(r)"), desc, fname="K[inhom]")

    # identify all close pairs
    rmax <- max(r)
    close <- closepairs(X, rmax)
    dIJ <- close$d
    # compute weights for these pairs
    I <- close$i
    J <- close$j
    wI <- reciplambda[I]
    wIJ <- 
      if(is.null(lambda2))
        reciplambda[I] * reciplambda[J]
      else 
        reciplambda2[cbind(I,J)]
    # 
    XI <- X[I]

    # compute edge corrected estimates
    if(any(correction == "border" | correction == "bord.modif")) {
      # border method
      # Compute distances to boundary
      b <- bdist.points(X)
      I <- close$i
      bI <- b[I]
      # apply reduced sample algorithm
      RS <- Kwtsum(dIJ, bI, wIJ, b, w=reciplambda, breaks)
      if(any(correction == "border")) {
        Kb <- RS$ratio
        if(renormalise) Kb <- Kb * renorm.factor
        K <- bind.fv(K, data.frame(border=Kb), "hat(%s^{bord})(r)",
                     "border-corrected estimate of %s",
                     "border")
      }
      if(any(correction == "bord.modif")) {
        Kbm <- RS$numerator/eroded.areas(W, r)
        if(renormalise) Kbm <- Kbm * renorm.factor
        K <- bind.fv(K, data.frame(bord.modif=Kbm), "hat(%s^{bordm})(r)",
                     "modified border-corrected estimate of %s",
                     "bord.modif")
      }
    }
    if(any(correction == "translate")) {
      # translation correction
      XJ <- X[J]
      edgewt <- edge.Trans(XI, XJ, paired=TRUE)
      allweight <- edgewt * wIJ
      wh <- whist(dIJ, breaks$val, allweight)
      Ktrans <- cumsum(wh)/area
      if(renormalise) Ktrans <- Ktrans * renorm.factor
      rmax <- diameter(W)/2
      Ktrans[r >= rmax] <- NA
      K <- bind.fv(K, data.frame(trans=Ktrans), "hat(%s^{trans})(r)",
                   "translation-correction estimate of %s",
                   "trans")
    }
    if(any(correction == "isotropic" | correction == "Ripley")) {
      # Ripley isotropic correction
      edgewt <- edge.Ripley(XI, matrix(dIJ, ncol=1))
      allweight <- edgewt * wIJ
      wh <- whist(dIJ, breaks$val, allweight)
      Kiso <- cumsum(wh)/area
      if(renormalise) Kiso <- Kiso * renorm.factor
      rmax <- diameter(W)/2
      Kiso[r >= rmax] <- NA
      K <- bind.fv(K, data.frame(iso=Kiso), "hat(%s^{iso})(r)",
                   "Ripley isotropic correction estimate of %s",
                   "iso")
    }

    # default is to display them all
    formula(K) <- . ~ r
    unitname(K) <- unitname(X)
    return(K)
}


Kwtsum <- function(dIJ, bI, wIJ, b, w, breaks) {
  #
  # "internal" routine to compute border-correction estimates of Kinhom
  #
  # dIJ:  vector containing pairwise distances for selected I,J pairs
  # bI:   corresponding vector of boundary distances for I
  # wIJ:  product weight for selected I, J pairs
  #
  # b:    vector of ALL distances to window boundary
  # w:   weights for ALL points
  #
  # breaks : breakpts object
  #

  stopifnot(length(dIJ) == length(bI))
  stopifnot(length(bI) == length(wIJ))
  stopifnot(length(w) == length(b))

  if(!is.finite(sum(w, wIJ)))
    stop("Weights in K-function were infinite or NA")
  
  # determine which distances d_{ij} were observed without censoring
  uncen <- (dIJ <= bI)
  #
  # histogram of noncensored distances
  nco <- whist(dIJ[uncen], breaks$val, wIJ[uncen])
  # histogram of censoring times for noncensored distances
  ncc <- whist(bI[uncen], breaks$val, wIJ[uncen])
  # histogram of censoring times (yes, this is a different total size)
  cen <- whist(b, breaks$val, w)
  # total weight of censoring times beyond rightmost breakpoint
  uppercen <- sum(w[b > max(breaks$val)])
  # go
  RS <- reduced.sample(nco, cen, ncc, show=TRUE, uppercen=uppercen)
  # extract results
  numerator   <- RS$numerator
  denominator <- RS$denominator
  ratio        <- RS$numerator/RS$denominator
  # check
  if(length(numerator) != breaks$ncells)
    stop("internal error: length(numerator) != breaks$ncells")
  if(length(denominator) != breaks$ncells)
    stop("internal error: length(denom.count) != breaks$ncells")
  return(list(numerator=numerator, denominator=denominator, ratio=ratio))
}
	
