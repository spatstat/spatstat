#
#	Kinhom.S	Estimation of K function for inhomogeneous patterns
#
#	$Revision: 1.72 $	$Date: 2014/02/11 08:32:33 $
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
  L <- eval.fv(sqrt(pmax.int(K,0)/pi))
  # relabel the fv object
  L <- rebadge.fv(L, quote(L[inhom](r)), c("L", "inhom"),
                  names(K), new.labl=attr(K, "labl"))
  attr(L, "labl") <- attr(K, "labl")
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
    diamW <- diameter(W)
    
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
                               good="good",
                               best="best"),
                             multi=TRUE)

    best.wanted <- ("best" %in% correction)
    ## replace 'good' by the optimal choice for this size of dataset
    if("good" %in% correction)
      correction[correction == "good"] <- good.correction.K(X)
    ## retain only corrections that are implemented for the window
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
  # Efficient code for border correction and no correction
  # Usable only if r values are evenly spaced from 0 to rmax
  # Invoked automatically if number of points is large

    can.do.fast <- breaks$even  && missing(lambda2)
    large.n    <- (npts >= nlarge)
    demand.best <- correction.given && best.wanted
    large.n.trigger <- large.n && !correction.given
    fastcorrections <- c("border", "bord.modif", "none")
    fastdefault <- "border"
    correction.fast  <- all(correction %in% fastcorrections)
    will.do.fast <- can.do.fast && (correction.fast || large.n.trigger)
    asked.fast <- (correction.given && correction.fast) ||
                  (nlarge.given && large.n.trigger)
    if(!can.do.fast && asked.fast) {
      whynot <-
        if(!(breaks$even)) "r values not evenly spaced" else
        if(!missing(lambda)) "matrix lambda2 was given" else NULL
      warning(paste(c("cannot use efficient code", whynot), sep="; "))
    }
    if(will.do.fast) {
      ## Compute Kinhom using fast algorithm(s)
      ## determine correction(s)
      ok <- correction %in% fastcorrections
      correction <- if(any(ok)) correction[ok] else fastdefault
      bord <- any(correction %in% c("border", "bord.modif"))
      none <- any(correction =="none")
      if(!all(ok)) {
        ## some corrections were overridden; notify user
        corx <- c(if(bord) "border correction estimate" else NULL,
                  if(none) "uncorrected estimate" else NULL)
        corx <- paste(corx, collapse=" and ")
        message(paste("number of data points exceeds",
                      nlarge, "- computing", corx , "only"))
      }
      ## restrict r values to recommended range, unless specifically requested
      if(!rfixed) 
        r <- seq(from=0, to=alim[2], length.out=length(r))
      ## border method
      if(bord) {
        Kb <- Kborder.engine(X, max(r), length(r), correction,
                             weights=reciplambda)
        Kb <- tweak.fv.entry(Kb, "border", new.labl="{hat(%s)[%s]^{bord}} (r)")
        Kb <- tweak.fv.entry(Kb, "bord.modif", new.labl="{hat(%s)[%s]^{bordm}} (r)")
      }
      ## uncorrected
      if(none) {
        Kn <- Knone.engine(X, max(r), length(r), weights=reciplambda)
        Kn <- tweak.fv.entry(Kn, "un", new.labl="{hat(%s)[%s]^{un}} (r)")
      }
      K <-
        if(bord && !none) Kb else
        if(!bord && none) Kn else 
      cbind.fv(Kb, Kn[, names(Kn) != "theo"])
      ## tweak labels
      K <- rebadge.fv(K, quote(K[inhom](r)), c("K", "inhom"))
      return(K)
    }

  ###########################################
  # Fast code for rectangular window
  ###########################################

  if(can.do.fast && is.rectangle(W) && spatstat.options("use.Krect")) {
    K <-  Krect.engine(X, rmax, length(r), correction,
                        weights=reciplambda, fname=c("K", "inhom"))
    K <- rebadge.fv(K, quote(K[inhom](r)), c("K", "inhom"))
    attr(K, "alim") <- alim
    return(K)
  }
  
  ###########################################
  # Slower code
  ###########################################
        
        
    # this will be the output data frame
    K <- data.frame(r=r, theo= pi * r^2)
    desc <- c("distance argument r", "theoretical Poisson %s")
    K <- fv(K, "r", quote(K[inhom](r)),
            "theo", , alim, c("r","{%s[%s]^{pois}}(r)"), desc,
            fname=c("K", "inhom"))

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
        K <- bind.fv(K, data.frame(border=Kb), "{hat(%s)[%s]^{bord}}(r)",
                     "border-corrected estimate of %s",
                     "border")
      }
      if(any(correction == "bord.modif")) {
        Kbm <- RS$numerator/eroded.areas(W, r)
        if(renormalise) Kbm <- Kbm * renorm.factor
        K <- bind.fv(K, data.frame(bord.modif=Kbm), "{hat(%s)[%s]^{bordm}}(r)",
                     "modified border-corrected estimate of %s",
                     "bord.modif")
      }
    }
    if(any(correction == "translate")) {
      # translation correction
      edgewt <- edge.Trans(dx=close$dx, dy=close$dy, W=W, paired=TRUE)
      allweight <- edgewt * wIJ
      wh <- whist(dIJ, breaks$val, allweight)
      Ktrans <- cumsum(wh)/area
      if(renormalise) Ktrans <- Ktrans * renorm.factor
      rmax <- diamW/2
      Ktrans[r >= rmax] <- NA
      K <- bind.fv(K, data.frame(trans=Ktrans), "{hat(%s)[%s]^{trans}}(r)",
                   "translation-correction estimate of %s",
                   "trans")
    }
    if(any(correction == "isotropic" | correction == "Ripley")) {
      # Ripley isotropic correction
      edgewt <- edge.Ripley(X[I], matrix(dIJ, ncol=1))
      allweight <- edgewt * wIJ
      wh <- whist(dIJ, breaks$val, allweight)
      Kiso <- cumsum(wh)/area
      if(renormalise) Kiso <- Kiso * renorm.factor
      rmax <- diamW/2
      Kiso[r >= rmax] <- NA
      K <- bind.fv(K, data.frame(iso=Kiso), "{hat(%s)[%s]^{iso}}(r)",
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
	
