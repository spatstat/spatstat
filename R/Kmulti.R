#
#	Kmulti.S		
#
#	Compute estimates of cross-type K functions
#	for multitype point patterns
#
#	$Revision: 5.49 $	$Date: 2018/10/13 06:55:41 $
#
#
# -------- functions ----------------------------------------
#	Kcross()	cross-type K function K_{ij}
#                       between types i and j
#
#	Kdot()          K_{i\bullet}
#                       between type i and all points regardless of type
#
#       Kmulti()        (generic)
#
#
# -------- standard arguments ------------------------------	
#	X		point pattern (of class 'ppp')
#				including 'marks' vector
#	r		distance values at which to compute K	
#
# -------- standard output ------------------------------
#      A data frame with columns named
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

"Lcross" <- function(X, i, j, ..., from, to) {
  if(!is.multitype(X, dfok=FALSE)) 
	stop("Point pattern must be multitype")
  if(missing(i)) i <- if(!missing(from)) from else levels(marks(X))[1]
  if(missing(j)) j <- if(!missing(to)) to else levels(marks(X))[2]
  K <- Kcross(X, i, j, ...)
  L <- eval.fv(sqrt(K/pi))
  # relabel the fv object
  iname <- make.parseable(paste(i))
  jname <- make.parseable(paste(j))
  L <- rebadge.fv(L,
                  substitute(L[i,j](r),
                             list(i=iname,j=jname)),
                  c("L", paste0("list(", iname, ",", jname, ")")),
                  new.yexp=substitute(L[list(i,j)](r),
                                      list(i=iname,j=jname)))
  attr(L, "labl") <- attr(K, "labl")
  return(L)  
}

"Ldot" <- function(X, i, ..., from) {
  if(!is.multitype(X, dfok=FALSE)) 
	stop("Point pattern must be multitype")
  if(missing(i)) i <- if(!missing(from)) from else levels(marks(X))[1]
  K <- Kdot(X, i, ...)
  L <- eval.fv(sqrt(K/pi))
  # relabel the fv object
  iname <- make.parseable(paste(i))
  L <- rebadge.fv(L,
                  substitute(L[i ~ dot](r), list(i=iname)),
                  c("L", paste(iname, "~ symbol(\"\\267\")")), 
                  new.yexp=substitute(L[i ~ symbol("\267")](r), list(i=iname)))
  attr(L, "labl") <- attr(K, "labl")
  return(L)  
}

"Kcross" <- 
function(X, i, j, r=NULL, breaks=NULL,
         correction =c("border", "isotropic", "Ripley", "translate") , ...,
         ratio=FALSE, from, to)
{
  verifyclass(X, "ppp")
  if(!is.multitype(X, dfok=FALSE)) 
	stop("Point pattern must be multitype")
  if(missing(correction))
    correction <- NULL
  marx <- marks(X)
  if(missing(i))
    i <- if(!missing(from)) from else levels(marx)[1]
  if(missing(j))
    j <- if(!missing(to)) to else levels(marx)[2]
  I <- (marx == i)
  if(!any(I))
    stop(paste("No points have mark i =", i))

  if(i == j) {
    ## use Kest
    result <- do.call(Kest,
                      resolve.defaults(list(X=X[I],
                                            r=r, breaks=breaks,
                                            correction=correction, ratio=ratio),
                                       list(rmax=NULL), ## forbidden 
                                       list(...)))
  } else {
    J <- (marx == j)
    if(!any(J))
      stop(paste("No points have mark j =", j))
    result <- Kmulti(X, I, J,
                     r=r, breaks=breaks,
                     correction=correction, ratio=ratio, ...)
  }
  iname <- make.parseable(paste(i))
  jname <- make.parseable(paste(j))
  result <-
    rebadge.fv(result, 
               substitute(Kcross[i,j](r), list(i=iname,j=jname)),
               c("K", paste0("list(", iname, ",", jname, ")")), 
               new.yexp=substitute(K[list(i,j)](r),
                                   list(i=iname,j=jname)))
  return(result)
}

"Kdot" <- 
function(X, i, r=NULL, breaks=NULL,
         correction = c("border", "isotropic", "Ripley", "translate") , ...,
         ratio=FALSE, from)
{
  verifyclass(X, "ppp")
  if(!is.multitype(X, dfok=FALSE)) 
	stop("Point pattern must be multitype")
  if(missing(correction))
    correction <- NULL

  marx <- marks(X)
  if(missing(i))
    i <- if(!missing(from)) from else levels(marx)[1]
        
  I <- (marx == i)
  J <- rep.int(TRUE, X$n)  # i.e. all points
	
  if(!any(I)) stop(paste("No points have mark i =", i))
	
  result <- Kmulti(X, I, J,
                   r=r, breaks=breaks, correction=correction, ..., ratio=ratio)
  iname <- make.parseable(paste(i))
  result <-
    rebadge.fv(result,
               substitute(K[i ~ dot](r), list(i=iname)),
               c("K", paste0(iname, "~ symbol(\"\\267\")")),
               new.yexp=substitute(K[i ~ symbol("\267")](r), list(i=iname)))
  return(result)
}


"Kmulti"<-
function(X, I, J, r=NULL, breaks=NULL,
         correction = c("border", "isotropic", "Ripley", "translate") , ...,
         ratio=FALSE)
{
  verifyclass(X, "ppp")

  npts <- npoints(X)
  W <- X$window
  areaW <- area(W)

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
                             best="best"),
                           multi=TRUE)

  correction <- implemented.for.K(correction, W$type, correction.given)

  I <- ppsubset(X, I)
  J <- ppsubset(X, J)
  if(is.null(I) || is.null(J))
    stop("I and J must be valid subset indices")
	
  if(!any(I)) stop("no points belong to subset I")
  if(!any(J)) stop("no points belong to subset J")
		
  nI <- sum(I)
  nJ <- sum(J)
  lambdaI <- nI/areaW
  lambdaJ <- nJ/areaW

  # r values 
  rmaxdefault <- rmax.rule("K", W, lambdaJ)
  breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max
        
  # recommended range of r values
  alim <- c(0, min(rmax, rmaxdefault))
        
  # this will be the output data frame
  # It will be given more columns later
  K <- data.frame(r=r, theo= pi * r^2)
  desc <- c("distance argument r", "theoretical Poisson %s")
  K <- fv(K, "r", quote(K[IJ](r)), 
          "theo", , alim, c("r","{%s[%s]^{pois}}(r)"),
          desc, fname=c("K", "list(I,J)"),
          yexp=quote(K[list(I,J)](r)))
  
  # save numerator and denominator?
  if(ratio) {
    denom <- lambdaI * lambdaJ * areaW
    numK <- eval.fv(denom * K)
    denK <- eval.fv(denom + K * 0)
    attributes(numK) <- attributes(denK) <- attributes(K)
    attr(numK, "desc")[2] <- "numerator for theoretical Poisson %s"
    attr(denK, "desc")[2] <- "denominator for theoretical Poisson %s"
  }

  # find close pairs of points
  XI <- X[I]
  XJ <- X[J]
  close <- crosspairs(XI, XJ, max(r), what="ijd")
# close$i and close$j are serial numbers in XI and XJ respectively;        
# map them to original serial numbers in X
  orig <- seq_len(npts)
  imap <- orig[I]
  jmap <- orig[J]
  iX <- imap[close$i]
  jX <- jmap[close$j]
# eliminate any identical pairs
  if(any(I & J)) {
    ok <- (iX != jX)
    if(!all(ok)) {
      close$i  <- close$i[ok]
      close$j  <- close$j[ok]
      close$d  <- close$d[ok]
    }
  }
# extract information for these pairs (relative to orderings of XI, XJ)
  dcloseIJ <- close$d
  icloseI  <- close$i
  jcloseJ  <- close$j
        
# Compute estimates by each of the selected edge corrections.
        
  if(any(correction == "none")) {
    # uncorrected! 
    wh <- whist(dcloseIJ, breaks$val)  # no weights
    numKun <- cumsum(wh)
    denKun <- lambdaI * lambdaJ * areaW
    Kun <- numKun/denKun
    K <- bind.fv(K, data.frame(un=Kun), "{hat(%s)[%s]^{un}}(r)",
                 "uncorrected estimate of %s",
                 "un")
    if(ratio) {
      # save numerator and denominator
      numK <- bind.fv(numK, data.frame(un=numKun), "{hat(%s)[%s]^{un}}(r)",
                 "numerator of uncorrected estimate of %s",
                 "un")
      denK <- bind.fv(denK, data.frame(un=denKun), "{hat(%s)[%s]^{un}}(r)",
                 "denominator of uncorrected estimate of %s",
                 "un")
    }

  }
  if(any(correction == "border" | correction == "bord.modif")) {
    # border method
    # distance to boundary from each point of type I
    bI <- bdist.points(XI)
    # distance to boundary from first element of each (i, j) pair
    bcloseI <- bI[icloseI]
    # apply reduced sample algorithm
    RS <- Kount(dcloseIJ, bcloseI, bI, breaks)
    if(any(correction == "bord.modif")) {
      denom.area <- eroded.areas(W, r)
      numKbm <- RS$numerator
      denKbm <- denom.area * nI * nJ
      Kbm <- numKbm/denKbm
      K <- bind.fv(K, data.frame(bord.modif=Kbm), "{hat(%s)[%s]^{bordm}}(r)",
                   "modified border-corrected estimate of %s",
                   "bord.modif")
      if(ratio) {
        # save numerator and denominator
        numK <- bind.fv(numK, data.frame(bord.modif=numKbm),
                        "{hat(%s)[%s]^{bordm}}(r)",
                        "numerator of modified border-corrected estimate of %s",
                        "bord.modif")
        denK <- bind.fv(denK, data.frame(bord.modif=denKbm),
                        "{hat(%s)[%s]^{bordm}}(r)",
                        "denominator of modified border-corrected estimate of %s",
                        "bord.modif")
      }
    }
    if(any(correction == "border")) {
      numKb <- RS$numerator
      denKb <- lambdaJ * RS$denom.count
      Kb <- numKb/denKb
      K <- bind.fv(K, data.frame(border=Kb), "{hat(%s)[%s]^{bord}}(r)",
                   "border-corrected estimate of %s",
                   "border")
      if(ratio) {
        numK <- bind.fv(numK, data.frame(border=numKb),
                        "{hat(%s)[%s]^{bord}}(r)",
                        "numerator of border-corrected estimate of %s",
                        "border")
        denK <- bind.fv(denK, data.frame(border=denKb),
                        "{hat(%s)[%s]^{bord}}(r)",
                        "denominator of border-corrected estimate of %s",
                        "border")
      }
    }
  }
  if(any(correction == "translate")) {
    # translation correction
    edgewt <- edge.Trans(XI[icloseI], XJ[jcloseJ], paired=TRUE)
    wh <- whist(dcloseIJ, breaks$val, edgewt)
    numKtrans <- cumsum(wh)
    denKtrans <- lambdaI * lambdaJ * areaW
    Ktrans <- numKtrans/denKtrans
    rmax <- diameter(W)/2
    Ktrans[r >= rmax] <- NA
    K <- bind.fv(K, data.frame(trans=Ktrans), "{hat(%s)[%s]^{trans}}(r)", 
                 "translation-corrected estimate of %s",
                 "trans")
    if(ratio) {
      numK <- bind.fv(numK, data.frame(trans=numKtrans),
                      "{hat(%s)[%s]^{trans}}(r)",
                      "numerator of translation-corrected estimate of %s",
                      "trans")
      denK <- bind.fv(denK, data.frame(trans=denKtrans),
                      "{hat(%s)[%s]^{trans}}(r)",
                      "denominator of translation-corrected estimate of %s",
                      "trans")
    }
  }
  if(any(correction == "isotropic")) {
    # Ripley isotropic correction
    edgewt <- edge.Ripley(XI[icloseI], matrix(dcloseIJ, ncol=1))
    wh <- whist(dcloseIJ, breaks$val, edgewt)
    numKiso <- cumsum(wh)
    denKiso <- lambdaI * lambdaJ * areaW
    Kiso <- numKiso/denKiso
    rmax <- diameter(W)/2
    Kiso[r >= rmax] <- NA
    K <- bind.fv(K, data.frame(iso=Kiso), "{hat(%s)[%s]^{iso}}(r)",
                 "Ripley isotropic correction estimate of %s",
                 "iso")
   if(ratio) {
      numK <- bind.fv(numK, data.frame(iso=numKiso), "{hat(%s)[%s]^{iso}}(r)",
                      "numerator of Ripley isotropic correction estimate of %s",
                      "iso")
      denK <- bind.fv(denK, data.frame(iso=denKiso), "{hat(%s)[%s]^{iso}}(r)",
                      "denominator of Ripley isotropic correction estimate of %s",
                      "iso")
    }
  }
  # default is to display them all
  formula(K) <- . ~ r
  unitname(K) <- unitname(X)
  
  if(ratio) {
    # finish up numerator & denominator
    formula(numK) <- formula(denK) <- . ~ r
    unitname(numK) <- unitname(denK) <- unitname(K)
    # tack on to result
    K <- rat(K, numK, denK, check=FALSE)
  }
  return(K)
}
