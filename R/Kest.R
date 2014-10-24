#
#	Kest.R		Estimation of K function
#
#	$Revision: 5.104 $	$Date: 2014/10/24 00:22:30 $
#
#
# -------- functions ----------------------------------------
#	Kest()		compute estimate of K
#                       using various edge corrections
#
#
# -------- standard arguments ------------------------------	
#	X		point pattern (of class 'ppp')
#
#	r		distance values at which to compute K	
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
#			using standard formula (denominator = count of points)
#
#       bord.modif:	K function estimated by border method
#			using modified formula 
#			(denominator = area of eroded window
#
# ------------------------------------------------------------------------

"Lest" <- function(X, ...) {
  K <- Kest(X, ...)
  L <- eval.fv(sqrt(K/pi))
  # handle variance estimates
  if(any(varcols <- colnames(K) %in% c("rip", "ls"))) {
    r <- with(L, .x)
    L[,varcols] <- as.data.frame(K)[,varcols]/(2 * pi * r)^2
    # fix 0/0
    n <- npoints(X)
    A <- area(Window(X))
    if(any(colnames(K) == "rip"))
      L[r == 0, "rip"] <- (2 * A/(n-1)^2)/(4 * pi)
    if(any(colnames(K) == "ls"))
      L[r == 0, "ls"]  <- (2 * A/(n * (n-1)))/(4 * pi)
  }
  # relabel the fv object
  L <- rebadge.fv(L, quote(L(r)), "L", names(K), new.labl=attr(K, "labl"))
  #
  return(L)  
}

"Kest"<-
function(X, ..., r=NULL, breaks=NULL, 
         correction=c("border", "isotropic", "Ripley", "translate"),
         nlarge=3000, domain=NULL, var.approx=FALSE,
         ratio=FALSE)
{
  verifyclass(X, "ppp")
  nlarge.given <- !missing(nlarge) && !is.null(nlarge)
  rfixed <- !is.null(r) || !is.null(breaks)
  npts <- npoints(X)
  W <- X$window
  areaW <- area(W)
  lambda <- npts/areaW
  lambda2 <- (npts * (npts - 1))/(areaW^2)

  if(!is.null(domain)) {
    # estimate based on contributions from a subdomain
    domain <- as.owin(domain)
    if(!is.subset.owin(domain, W))
      stop(paste(dQuote("domain"),
                 "is not a subset of the window of X"))
    # trick Kdot() into doing it
    indom <- factor(inside.owin(X$x, X$y, domain), levels=c(FALSE,TRUE))
    Kd <- Kdot(X %mark% indom, i="TRUE",
               r=r, breaks=breaks, correction=correction,
               ratio=ratio)
    # relabel and exit
    Kd <- rebadge.fv(Kd, quote(K(r)), "K")
    return(Kd)
  }

  rmaxdefault <- rmax.rule("K", W, lambda)        
  breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max

  # choose correction(s)
  correction.given <- !missing(correction) && !is.null(correction)
  if(is.null(correction))
    correction <- c("border", "isotropic", "Ripley", "translate")
  correction <- pickoption("correction", correction,
                           c(none="none",
                             border="border",
                             "bord.modif"="bord.modif",
                             isotropic="isotropic",
                             Ripley="isotropic",
                             trans="translate",
                             translate="translate",
                             translation="translate",
                             rigid="rigid",
                             good="good",
                             best="best"),
                           multi=TRUE)
  best.wanted <- ("best" %in% correction)
  # replace 'good' by the optimal choice for this size of dataset
  if("good" %in% correction)
    correction[correction == "good"] <- good.correction.K(X)
  # retain only corrections that are implemented for the window
  correction <- implemented.for.K(correction, W$type, correction.given)
  
  # recommended range of r values
  alim <- c(0, min(rmax, rmaxdefault))

  ###########################################
  # Efficient code for border correction and no correction
  # Usable only if r values are evenly spaced from 0 to rmax
  # Invoked automatically if number of points is large

  can.do.fast <- breaks$even
  large.n    <- (npts >= nlarge)
  demand.best <- correction.given && best.wanted
  large.n.trigger <- large.n && !correction.given
  fastcorrections <- c("border", "bord.modif", "none")
  fastdefault     <- "border"
  correction.fast   <- all(correction %in% fastcorrections)
  will.do.fast <- can.do.fast && (correction.fast || large.n.trigger)
  asked <- correction.fast || (nlarge.given && large.n.trigger)
  if(asked && !can.do.fast)
    warning("r values not evenly spaced - cannot use efficient code")
  if(will.do.fast) {
    # determine correction(s)
    ok <- correction %in% fastcorrections
    correction <- if(any(ok)) correction[ok] else fastdefault
    bord <- any(correction %in% c("border", "bord.modif"))
    none <- any(correction =="none")
    if(!all(ok)) {
      # some corrections were overridden; notify user
      corx <- c(if(bord) "border correction estimate" else NULL,
                if(none) "uncorrected estimate" else NULL)
      corx <- paste(corx, collapse=" and ")
      message(paste("number of data points exceeds",
                    nlarge, "- computing", corx , "only"))
    }
    # restrict r values to recommended range, unless specifically requested
    if(!rfixed) 
      r <- seq(from=0, to=alim[2], length.out=length(r))
    if(bord)
      Kb <- Kborder.engine(X, max(r), length(r), correction, ratio=ratio)
    if(none)
      Kn <- Knone.engine(X, max(r), length(r), ratio=ratio)
    if(bord && none) 
      return(cbind.fv(Kb, Kn[, names(Kn) != "theo"]))
    if(bord) return(Kb)
    if(none) return(Kn) 
  }

  do.fast.rectangle <-
    can.do.fast && is.rectangle(W) &&
      spatstat.options("use.Krect") && !any(correction == "rigid")
  
  if(do.fast.rectangle) {
    ###########################################
    ## Fast code for rectangular window
    ###########################################
    K <-  Krect.engine(X, rmax, length(r), correction, ratio=ratio)
    attr(K, "alim") <- alim
  } else {
    ###########################################
    ## Slower code
    ###########################################

    ## this will be the output data frame
    Kdf <- data.frame(r=r, theo = pi * r^2)
    desc <- c("distance argument r", "theoretical Poisson %s")
    denom <- lambda2 * areaW
    K <- ratfv(Kdf, NULL, denom,
               "r", quote(K(r)),
               "theo", NULL, alim, c("r","%s[pois](r)"), desc, fname="K",
               ratio=ratio)
  
    ## identify all close pairs
    rmax <- max(r)
    close <- closepairs(X, rmax)
    DIJ <- close$d

    if(any(correction == "none")) {
      ## uncorrected! For demonstration purposes only!
      wh <- whist(DIJ, breaks$val)  # no weights
      numKun <- cumsum(wh)
      denKun <- lambda2 * areaW
      ## uncorrected estimate of K
      K <- bind.ratfv(K,
                      data.frame(un=numKun), denKun,
                      "hat(%s)[un](r)",
                      "uncorrected estimate of %s",
                      "un",
                      ratio=ratio)
    }
  
    if(any(correction == "border" | correction == "bord.modif")) {
      ## border method
      ## Compute distances to boundary
      b <- bdist.points(X)
      I <- close$i
      bI <- b[I]
      ## apply reduced sample algorithm
      RS <- Kount(DIJ, bI, b, breaks)
      if(any(correction == "bord.modif")) {
        ## modified border correction
        denom.area <- eroded.areas(W, r)
        numKbm <- RS$numerator
        denKbm <- lambda2 * denom.area
        K <- bind.ratfv(K,
                        data.frame(bord.modif=numKbm),
                        data.frame(bord.modif=denKbm),
                        "hat(%s)[bordm](r)",
                        "modified border-corrected estimate of %s",
                        "bord.modif",
                        ratio=ratio)
      }
      if(any(correction == "border")) {
        numKb <- RS$numerator
        denKb <- lambda * RS$denom.count
        K <- bind.ratfv(K,
                        data.frame(border=numKb), 
                        data.frame(border=denKb), 
                        "hat(%s)[bord](r)",
                        "border-corrected estimate of %s",
                        "border",
                        ratio=ratio)
      }
    }

    if(any(correction == "translate")) {
      ## Ohser-Stoyan translation correction
      edgewt <- edge.Trans(dx=close$dx, dy=close$dy, W=W, paired=TRUE)
      wh <- whist(DIJ, breaks$val, edgewt)
      numKtrans <- cumsum(wh)
      denKtrans <- lambda2 * areaW
      h <- diameter(as.rectangle(W))/2
      numKtrans[r >= h] <- NA
      K <- bind.ratfv(K,
                      data.frame(trans=numKtrans),
                      denKtrans,
                      "hat(%s)[trans](r)",
                      "translation-corrected estimate of %s",
                      "trans",
                      ratio=ratio)
    }
    if(any(correction == "rigid")) {
      ## Ohser-Stoyan rigid motion correction
      CW <- rotmean(setcov(W))
      edgewt <- areaW/as.function(CW)(DIJ)
      wh <- whist(DIJ, breaks$val, edgewt)
      numKrigid <- cumsum(wh)
      denKrigid <- lambda2 * areaW
      h <- diameter(as.rectangle(W))
      numKrigid[r >= h] <- NA
      K <- bind.ratfv(K,
                      data.frame(rigid=numKrigid),
                      denKrigid,
                      "hat(%s)[rigid](r)",
                      "rigid motion-corrected estimate of %s",
                      "rigid",
                      ratio=ratio)
    }
    if(any(correction == "isotropic")) {
      ## Ripley isotropic correction
      XI <- ppp(close$xi, close$yi, window=W, check=FALSE)
      edgewt <- edge.Ripley(XI, matrix(DIJ, ncol=1))
      wh <- whist(DIJ, breaks$val, edgewt)
      numKiso <- cumsum(wh)
      denKiso <- lambda2 * areaW
      h <- diameter(W)/2
      numKiso[r >= h] <- NA
      K <- bind.ratfv(K,
                      data.frame(iso=numKiso),
                      denKiso,
                      "hat(%s)[iso](r)",
                      "Ripley isotropic correction estimate of %s",
                      "iso",
                      ratio=ratio)
    }
  }

  #############################
  ##  VARIANCE APPROXIMATION
  #############################

  if(var.approx) {
    ## Compute variance approximations
    A <- areaW
    P <- perimeter(W)
    n <- npts
    ## Ripley asymptotic approximation
    rip <- 2 * ((A/(n-1))^2) * (pi * r^2/A + 0.96 * P * r^3/A^2
                                + 0.13 * (n/A) * P * r^5/A^2)
    if(!ratio) {
      K <- bind.fv(K, data.frame(rip=rip),
                 "vR(r)", 
                 "Ripley approximation to var(%s) under CSR",
                 "iso")
    } else {
      den <- (n-1)^2
      ripnum <- den * rip
      ripden <- rep.int(den, length(rip))
      K <- bind.ratfv(K,
                      data.frame(rip=ripnum),
                      data.frame(rip=ripden),
                      "vR(r)", 
                      "Ripley approximation to var(%s) under CSR",
                      "iso")
    }
    if(W$type == "rectangle") {
      # Lotwick-Silverman
      a1r <- (0.21 * P * r^3 + 1.3 * r^4)/A^2
      a2r <- (0.24 * P * r^5 + 2.62 * r^6)/A^3
      # contains correction to typo on p52 of Diggle 2003
      # cf Lotwick & Silverman 1982 eq (5)
      br <- (pi * r^2/A) * (1 - pi * r^2/A) +
        (1.0716 * P * r^3 + 2.2375 * r^4)/A^2
      ls <- (A^2) * (2 * br - a1r + (n-2) * a2r)/(n*(n-1))
      # add column 
      if(!ratio) {
        K <- bind.fv(K, data.frame(ls=ls), "vLS(r)",
                     "Lotwick-Silverman approx to var(%s) under CSR",
                     "iso")
      } else {
        den <- n*(n-1)
        lsnum <- ls * den
        lsden <- rep.int(den, length(ls))
        K <- bind.ratfv(K,
                        data.frame(ls=lsnum),
                        data.frame(ls=lsden),
                        "vLS(r)",
                        "Lotwick-Silverman approx to var(%s) under CSR",
                        "iso")
      }
    }
  }

  ### FINISH OFF #####
  ## default plot will display all edge corrections
  formula(K) <- . ~ r
  nama <- rev(colnames(K))
  nama <- nama[!(nama %in% c("r", "rip", "ls"))]
  fvnames(K, ".") <- nama
  unitname(K) <- unitname(X)
  # copy to other components
  if(ratio)
    K <- conform.ratfv(K)

  return(K)
}

################################################################  
#############  SUPPORTING ALGORITHMS ###########################
################################################################  

Kount <- function(dIJ, bI, b, breaks) {
  #
  # "internal" routine to compute border-correction estimate of K or Kij
  #
  # dIJ:  vector containing pairwise distances for selected I,J pairs
  # bI:   corresponding vector of boundary distances for I
  # b:    vector of ALL distances to window boundary
  #
  # breaks : breakpts object
  #

  stopifnot(length(dIJ) == length(bI))
  
  # determine which distances d_{ij} were observed without censoring
  uncen <- (dIJ <= bI)
  # histogram of noncensored distances
  nco <- whist(dIJ[uncen], breaks$val)
  # histogram of censoring times for noncensored distances
  ncc <- whist(bI[uncen], breaks$val)
  # histogram of censoring times (yes, this is a different total size)
  cen <- whist(b, breaks$val)
  # count censoring times beyond rightmost breakpoint
  uppercen <- sum(b > max(breaks$val))
  # go
  RS <- reduced.sample(nco, cen, ncc, show=TRUE, uppercen=uppercen)
  # extract results
  numerator <- RS$numerator
  denom.count <- RS$denominator
  # check
  if(length(numerator) != breaks$ncells)
    stop("internal error: length(numerator) != breaks$ncells")
  if(length(denom.count) != breaks$ncells)
    stop("internal error: length(denom.count) != breaks$ncells")
  
  return(list(numerator=numerator, denom.count=denom.count))
}

#### interface to C code for border method

Kborder.engine <- function(X, rmax, nr=100,
                           correction=c("border", "bord.modif"),
                           weights=NULL, ratio=FALSE) 
{
  verifyclass(X, "ppp")
  npts <- npoints(X)
  W <- as.owin(X)

  areaW <- area(W)
  lambda <- npts/areaW
  lambda2 <- (npts * (npts - 1))/(areaW^2)

  if(missing(rmax))
    rmax <- diameter(W)/4
  r <- seq(from=0, to=rmax, length.out=nr)

  # this will be the output data frame
  Kdf <- data.frame(r=r, theo= pi * r^2)
  desc <- c("distance argument r", "theoretical Poisson %s")
  Kfv <- fv(Kdf, "r", quote(K(r)),
          "theo", , c(0,rmax), c("r","%s[pois](r)"), desc, fname="K")

  if(ratio) {
    # save numerator and denominator
    denom <- lambda2 * areaW
    numK <- eval.fv(denom * Kfv)
    denK <- eval.fv(denom + Kfv * 0)
    attributes(numK) <- attributes(denK) <- attributes(Kfv)
    numK <- rebadge.fv(numK, tags="theo",
                       new.desc="numerator for theoretical Poisson %s")
    denK <- rebadge.fv(denK, tags="theo",
                       new.desc="denominator for theoretical Poisson %s")
  }
  
  ####### start computing ############
  # sort in ascending order of x coordinate
  orderX <- fave.order(X$x)
  Xsort <- X[orderX]
  x <- Xsort$x
  y <- Xsort$y
  
  # boundary distances
  b <- bdist.points(Xsort)

  # call the C code
  if(is.null(weights)) {
    # determine whether the numerator can be stored as an integer
    bigint <- .Machine$integer.max
    if(npts < sqrt(bigint)) {
      # yes - use faster integer arithmetic
      res <- .C("KborderI",
                nxy=as.integer(npts),
                x=as.double(x),
                y=as.double(y),
                b=as.double(b),
                nr=as.integer(nr),
                rmax=as.double(rmax),
                numer=as.integer(integer(nr)),
                denom=as.integer(integer(nr)))
    } else {
      # no - need double precision storage
      res <- .C("KborderD",
                nxy=as.integer(npts),
                x=as.double(x),
                y=as.double(y),
                b=as.double(b),
                nr=as.integer(nr),
                rmax=as.double(rmax),
                numer=as.double(numeric(nr)),
                denom=as.double(numeric(nr)))
    }
    if("bord.modif" %in% correction) {
      denom.area <- eroded.areas(W, r)
      numKbm <- res$numer
      denKbm <- lambda2 * denom.area
      bm <- numKbm/denKbm
      Kfv <- bind.fv(Kfv, data.frame(bord.modif=bm), "hat(%s)[bordm](r)",
                   "modified border-corrected estimate of %s",
                   "bord.modif")
      if(ratio) {
        # save numerator and denominator
        numK <- bind.fv(numK, data.frame(bord.modif=numKbm),
                        "hat(%s)[bordm](r)",
                        "numerator of modified border-corrected estimate of %s",
                        "bord.modif")
        denK <- bind.fv(denK, data.frame(bord.modif=denKbm),
                        "hat(%s)[bordm](r)",
                        "denominator of modified border-corrected estimate of %s",
                        "bord.modif")
      }
    }
    if("border" %in% correction) {
      numKb <- res$numer
      denKb <- lambda * res$denom
      bord <- numKb/denKb
      Kfv <- bind.fv(Kfv, data.frame(border=bord), "hat(%s)[bord](r)",
                   "border-corrected estimate of %s",
                   "border")
      if(ratio) {
        numK <- bind.fv(numK, data.frame(border=numKb),
                        "hat(%s)[bord](r)",
                        "numerator of border-corrected estimate of %s",
                        "border")
        denK <- bind.fv(denK, data.frame(border=denKb),
                        "hat(%s)[bord](r)",
                        "denominator of border-corrected estimate of %s",
                        "border")
      }
    }
  } else {
    # weighted version
    if(is.numeric(weights)) {
      if(length(weights) != X$n)
        stop("length of weights argument does not match number of points in X")
    } else {
      wim <- as.im(weights, W)
      weights <- wim[X, drop=FALSE]
      if(any(is.na(weights)))
        stop("domain of weights image does not contain all points of X")
    }
    weights.Xsort <- weights[orderX]
    res <- .C("Kwborder",
              nxy=as.integer(npts),
              x=as.double(x),
              y=as.double(y),
              w=as.double(weights.Xsort),
              b=as.double(b),
              nr=as.integer(nr),
              rmax=as.double(rmax),
              numer=as.double(numeric(nr)),
              denom=as.double(numeric(nr)))
    if("border" %in% correction) {
      bord <- res$numer/res$denom
      Kfv <- bind.fv(Kfv, data.frame(border=bord), "hat(%s)[bord](r)",
                     "border-corrected estimate of %s",
                     "border")
      if(ratio) {
        numK <- bind.fv(numK, data.frame(border=res$numer),
                        "hat(%s)[bord](r)",
                        "numerator of border-corrected estimate of %s",
                        "border")
        denK <- bind.fv(denK, data.frame(border=res$denom),
                        "hat(%s)[bord](r)",
                        "denominator of border-corrected estimate of %s",
                        "border")
      }
    }
    if("bord.modif" %in% correction) {
      numKbm <- res$numer
      denKbm <- eroded.areas(W, r)
      bm <- numKbm/denKbm
      Kfv <- bind.fv(Kfv, data.frame(bord.modif=bm), "hat(%s)[bordm](r)",
                     "modified border-corrected estimate of %s",
                     "bord.modif")
      if(ratio) {
        # save numerator and denominator
        numK <- bind.fv(numK, data.frame(bord.modif=numKbm),
                        "hat(%s)[bordm](r)",
                        "numerator of modified border-corrected estimate of %s",
                        "bord.modif")
        denK <- bind.fv(denK, data.frame(bord.modif=denKbm),
                        "hat(%s)[bordm](r)",
                        "denominator of modified border-corrected estimate of %s",
                        "bord.modif")
      }
    }
  }
  ##
  # default is to display them all
  formula(Kfv) <- . ~ r
  unitname(Kfv) <- unitname(X)
  if(ratio) {
    # finish off numerator and denominator
    formula(numK) <- formula(denK) <- . ~ r
    unitname(denK) <- unitname(numK) <- unitname(X)
    # tack on to result
    Kfv <- rat(Kfv, numK, denK, check=FALSE)
  }
  return(Kfv)
}

Knone.engine <- function(X, rmax, nr=100,
                         weights=NULL, ratio=FALSE) 
{
  verifyclass(X, "ppp")
  npts <- npoints(X)
  W <- as.owin(X)

  areaW <- area(W)
  lambda <- npts/areaW
  lambda2 <- (npts * (npts - 1))/(areaW^2)
  denom <- lambda2 * areaW

  if(missing(rmax))
    rmax <- diameter(W)/4
  r <- seq(from=0, to=rmax, length.out=nr)

  # this will be the output data frame
  Kdf <- data.frame(r=r, theo= pi * r^2)
  desc <- c("distance argument r", "theoretical Poisson %s")
  Kfv <- fv(Kdf, "r", quote(K(r)),
          "theo", , c(0,rmax), c("r","%s[pois](r)"), desc, fname="K")

  if(ratio) {
    # save numerator and denominator
    numK <- eval.fv(denom * Kfv)
    denK <- eval.fv(denom + Kfv * 0)
    attributes(numK) <- attributes(denK) <- attributes(Kfv)
    numK <- rebadge.fv(numK, tags="theo",
                       new.desc="numerator for theoretical Poisson %s")
    denK <- rebadge.fv(denK, tags="theo",
                       new.desc="denominator for theoretical Poisson %s")
  }
  
  ####### start computing ############
  # sort in ascending order of x coordinate
  orderX <- fave.order(X$x)
  Xsort <- X[orderX]
  x <- Xsort$x
  y <- Xsort$y
  
  # call the C code
  if(is.null(weights)) {
    # determine whether the numerator can be stored as an integer
    bigint <- .Machine$integer.max
    if(npts < sqrt(bigint)) {
      # yes - use faster integer arithmetic
      res <- .C("KnoneI",
                nxy=as.integer(npts),
                x=as.double(x),
                y=as.double(y),
                nr=as.integer(nr),
                rmax=as.double(rmax),
                numer=as.integer(integer(nr)))
    } else {
      # no - need double precision storage
      res <- .C("KnoneD",
                nxy=as.integer(npts),
                x=as.double(x),
                y=as.double(y),
                nr=as.integer(nr),
                rmax=as.double(rmax),
                numer=as.double(numeric(nr)))
    }

    numKun <- res$numer
    denKun <- denom # = lambda2 * areaW
    Kun <- numKun/denKun
  } else {
    # weighted version
    if(is.numeric(weights)) {
      if(length(weights) != X$n)
        stop("length of weights argument does not match number of points in X")
    } else {
      wim <- as.im(weights, W)
      weights <- wim[X, drop=FALSE]
      if(any(is.na(weights)))
        stop("domain of weights image does not contain all points of X")
    }
    weights.Xsort <- weights[orderX]
    res <- .C("Kwnone",
              nxy=as.integer(npts),
              x=as.double(x),
              y=as.double(y),
              w=as.double(weights.Xsort),
              nr=as.integer(nr),
              rmax=as.double(rmax),
              numer=as.double(numeric(nr)))
    numKun <- res$numer
    denKun <- sum(weights)
    Kun <- numKun/denKun
  }

  # tack on to fv object
  Kfv <- bind.fv(Kfv, data.frame(un=Kun), "hat(%s)[un](r)",
                 "uncorrected estimate of %s",
                 "un")
  if(ratio) {
    numK <- bind.fv(numK, data.frame(un=numKun),
                    "hat(%s)[un](r)",
                    "numerator of uncorrected estimate of %s",
                    "un")
    denK <- bind.fv(denK, data.frame(un=denKun),
                    "hat(%s)[un](r)",
                    "denominator of uncorrected estimate of %s",
                    "un")
  }
  ##
  # default is to display them all
  formula(Kfv) <- . ~ r
  unitname(Kfv) <- unitname(X)
  if(ratio) {
    # finish off numerator and denominator
    formula(numK) <- formula(denK) <- . ~ r
    unitname(denK) <- unitname(numK) <- unitname(X)
    # tack on to result
    Kfv <- rat(Kfv, numK, denK, check=FALSE)
  }
  return(Kfv)
}

     

rmax.rule <- function(fun="K", W, lambda) {
  verifyclass(W, "owin")
  switch(fun,
         K = {
           # Ripley's Rule
           ripley <- min(diff(W$xrange), diff(W$yrange))/4
           # Count at most 1000 neighbours per point
           rlarge <- if(!missing(lambda)) sqrt(1000 /(pi * lambda)) else Inf
           rmax <- min(rlarge, ripley)
         },
         F = ,
         G = ,
         J = {
           # rule of thumb
           rdiam  <- diameter(as.rectangle(W))/2
           # Poisson process has F(rlarge) = 1 - 10^(-5)
           rlarge <-
             if(!missing(lambda)) sqrt(log(1e5)/(pi * lambda)) else Inf
           rmax <- min(rlarge, rdiam)
         },
         stop(paste("Unrecognised function type", sQuote(fun)))
         )
  return(rmax)
}
           
    
implemented.for.K <- function(correction, windowtype, explicit) {
  pixels <- (windowtype == "mask")
  if(any(correction == "best")) {
    # select best available correction
    correction <- if(!pixels) "isotropic" else "translate"
  } else {
    # available selection of edge corrections depends on window
    if(pixels) {
      iso <- (correction == "isotropic") 
      if(any(iso)) {
        whinge <- "Isotropic correction not implemented for binary masks"
        if(explicit) {
          if(all(iso)) stop(whinge) else warning(whinge)
        }
        correction <- correction[!iso]
      }
    }
  }
  return(correction)
}

good.correction.K <- function(X) {
  nX <- npoints(X)
  W <- as.owin(X)
  avail <- c("none",
             if(nX < 1e5) "border" else NULL,
             if(nX < 3000)"translate" else NULL,
             if(nX < 1000 && !is.mask(W)) "isotropic" else NULL)
  chosen <- rev(avail)[1]
  return(chosen)
}

Krect.engine <- function(X, rmax, nr=100,
                           correction,
                           weights=NULL, ratio=FALSE, fname="K") {
  verifyclass(X, "ppp")
  npts <- npoints(X)
  W <- as.owin(X)

  areaW <- area(W)
  width <- sidelengths(W)[1]
  height <- sidelengths(W)[2]
  lambda <- npts/areaW
  lambda2 <- (npts * (npts - 1))/(areaW^2)

  if(missing(rmax))
    rmax <- diameter(W)/4
  r <- seq(from=0, to=rmax, length.out=nr)

  if(weighted <- !is.null(weights)) {
    ## coerce weights to a vector
    if(is.numeric(weights)) {
      check.nvector(weights, npts)
    } else {
      wim <- as.im(weights, W)
      weights <- wim[X, drop=FALSE]
      if(any(is.na(weights)))
        stop("domain of weights image does not contain all points of X")
    }
  }

  # this will be the output data frame
  Kdf <- data.frame(r=r, theo= pi * r^2)
  desc <- c("distance argument r", "theoretical Poisson %s")
  denom <- if(weighted) areaW else (lambda2 * areaW)
  Kfv <- ratfv(Kdf, NULL, denom,
               "r", quote(K(r)),
               "theo", NULL, c(0,rmax),
               c("r", makefvlabel(NULL, NULL, fname, "pois")),
               desc, fname=fname,
               ratio=ratio)

  ####### prepare data ############

  if(!all(correction == "translate")) {
    ## Ensure rectangle has its bottom left corner at the origin
    if(W$xrange[1] != 0 || W$yrange[1] != 0) {
      X <- shift(X, origin="bottomleft")
      W <- as.owin(X)
    }
  }

  ## sort in ascending order of x coordinate
  orderX <- fave.order(X$x)
  x <- X$x[orderX]
  y <- X$y[orderX]
  if(weighted)
    wt <- weights[orderX]

  ## establish algorithm parameters
  doIso <- "isotropic" %in% correction 
  doTrans <- "translate" %in% correction
  doBord <- any(c("border", "bord.modif") %in% correction)
  doUnco <- "none" %in% correction
  trimedge <- spatstat.options("maxedgewt")

  ## allocate space for results
  ziso   <- numeric(if(doIso) nr else 1L)
  ztrans <- numeric(if(doTrans) nr else 1L)
  
  ## call the C code
  if(weighted) {
    ## weighted version
    zbnumer <- numeric(if(doBord) nr else 1L)
    zbdenom <- numeric(if(doBord) nr else 1L)
    zunco   <- numeric(if(doUnco) nr else 1L)
    res <- .C("KrectWtd",
              width=as.double(width),
              height=as.double(height),
              nxy=as.integer(npts),
              x=as.double(x),
              y=as.double(y),
              w=as.double(wt),
              nr=as.integer(nr),
              rmax=as.double(rmax),
              trimedge=as.double(trimedge),
              doIso=as.integer(doIso),
              doTrans=as.integer(doTrans),
              doBord=as.integer(doBord),
              doUnco=as.integer(doUnco),
              iso=as.double(ziso),
              trans=as.double(ztrans),
              bnumer=as.double(zbnumer),
              bdenom=as.double(zbdenom),
              unco=as.double(zunco))
  } else if(npts < sqrt(.Machine$integer.max)) {
    ## unweighted
    ## numerator of border correction can be stored as an integer
    ## use faster integer arithmetic
    zbnumer <- integer(if(doBord) nr else 1L)
    zbdenom <- integer(if(doBord) nr else 1L)
    zunco   <- integer(if(doUnco) nr else 1L)
    res <- .C("KrectInt",
              width=as.double(width),
              height=as.double(height),
              nxy=as.integer(npts),
              x=as.double(x),
              y=as.double(y),
              nr=as.integer(nr),
              rmax=as.double(rmax),
              trimedge=as.double(trimedge),
              doIso=as.integer(doIso),
              doTrans=as.integer(doTrans),
              doBord=as.integer(doBord),
              doUnco=as.integer(doUnco),
              iso=as.double(ziso),
              trans=as.double(ztrans),
              bnumer=as.integer(zbnumer),
              bdenom=as.integer(zbdenom),
              unco=as.integer(zunco))
  } else {
    ## unweighted
    ## need double precision storage
    zbnumer <- numeric(if(doBord) nr else 1L)
    zbdenom <- numeric(if(doBord) nr else 1L)
    zunco   <- numeric(if(doUnco) nr else 1L)
    res <- .C("KrectDbl",
              width=as.double(width),
              height=as.double(height),
              nxy=as.integer(npts),
              x=as.double(x),
              y=as.double(y),
              nr=as.integer(nr),
              rmax=as.double(rmax),
              trimedge=as.double(trimedge),
              doIso=as.integer(doIso),
              doTrans=as.integer(doTrans),
              doBord=as.integer(doBord),
              doUnco=as.integer(doUnco),
              iso=as.double(ziso),
              trans=as.double(ztrans),
              bnumer=as.double(zbnumer),
              bdenom=as.double(zbdenom),
              unco=as.double(zunco))
  }

  ## Process corrections in reverse order of priority

  ## Uncorrected estimate
  if("none" %in% correction) {
    numKun <- res$unco
    denKun <- if(weighted) areaW else (lambda2 * areaW)
    Kfv <- bind.ratfv(Kfv,
                      data.frame(un=numKun),
                      denKun,
                      makefvlabel(NULL, "hat", fname, "un"),
                      "uncorrected estimate of %s",
                      "un",
                      ratio=ratio)
  }
  
  ## Modified border correction
  if("bord.modif" %in% correction) {
    denom.area <- eroded.areas(W, r)
    numKbm <- res$bnumer
    denKbm <- if(weighted) denom.area else (lambda2 * denom.area)
    Kfv <- bind.ratfv(Kfv,
                      data.frame(bord.modif=numKbm),
                      denKbm,
                      makefvlabel(NULL, "hat", fname, "bordm"),
                      "modified border-corrected estimate of %s",
                      "bord.modif",
                      ratio=ratio)
  }
  ## Border correction
  if("border" %in% correction) {
    numKb <- res$bnumer
    denKb <- if(weighted) res$bdenom else lambda * res$bdenom
    Kfv <- bind.ratfv(Kfv,
                      data.frame(border=numKb),
                      denKb,
                      makefvlabel(NULL, "hat", fname, "bord"),
                      "border-corrected estimate of %s",
                      "border",
                      ratio=ratio)
  }
  
  ## translation correction
  if("translate" %in% correction) {
    numKtrans <- res$trans
    denKtrans <- if(weighted) areaW else (lambda2 * areaW)
    h <- diameter(as.rectangle(W))/2
    numKtrans[r >= h] <- NA
    Kfv <- bind.ratfv(Kfv,
                      data.frame(trans=numKtrans),
                      denKtrans,
                      makefvlabel(NULL, "hat", fname, "trans"),
                      "translation-corrected estimate of %s",
                      "trans",
                      ratio=ratio)
  }
  ## isotropic correction
  if("isotropic" %in% correction) {
    numKiso <- res$iso
    denKiso <- if(weighted) areaW else (lambda2 * areaW)
    h <- diameter(as.rectangle(W))/2
    numKiso[r >= h] <- NA
    Kfv <- bind.ratfv(Kfv,
                      data.frame(iso=numKiso),
                      denKiso,
                      makefvlabel(NULL, "hat", fname, "iso"),
                      "isotropic-corrected estimate of %s",
                      "iso",
                      ratio=ratio)
  }
  ##
  # default is to display them all
  formula(Kfv) <- . ~ r
  unitname(Kfv) <- unitname(X)
  if(ratio) 
    Kfv <- conform.ratfv(Kfv)
  return(Kfv)
}



