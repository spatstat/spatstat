#
#	Ksector.R	Estimation of 'sector K function'
#
#	$Revision: 1.1 $	$Date: 2014/07/29 12:04:50 $
#

Ksector <- function(X, begin=0, end=2*pi, ..., r=NULL, breaks=NULL, 
                    correction=c("border", "isotropic", "Ripley", "translate"),
                    ratio=FALSE)
{
  verifyclass(X, "ppp")
  rfixed <- !is.null(r) || !is.null(breaks)
  npts <- npoints(X)
  W <- X$window
  area <- area.owin(W)
  lambda <- npts/area
  lambda2 <- (npts * (npts - 1))/(area^2)
  rmaxdefault <- rmax.rule("K", W, lambda)        
  breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max

  check.1.real(begin)
  check.1.real(end)
  check.in.range(begin, c(0, 2*pi))
  check.in.range(end, c(0, 2*pi))
  stopifnot(begin < end)

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

  # this will be the output data frame
  Kdf <- data.frame(r=r, theo = ((end-begin)/2) * r^2)
  desc <- c("distance argument r", "theoretical Poisson %s")
  denom <- lambda2 * area
  K <- ratfv(Kdf, NULL, denom,
             "r", quote(K(r)),
             "theo", NULL, alim,
             c("r","{%s[%s]^{pois}}(r)"),
             desc,
             fname=c("K", "sector"),
             ratio=ratio)

  # identify all close pairs
  rmax <- max(r)
  close <- as.data.frame(closepairs(X, rmax))

  ## select pairs in angular range
  ang <- with(close, atan2(dy, dx)) %% (2*pi)
  ok <- (begin <= ang) & (ang <= end)
  close <- close[ok, , drop=FALSE]

  ## pairwise distances
  DIJ <- close$d

  if(any(correction == "none")) {
    # uncorrected! For demonstration purposes only!
    wh <- whist(DIJ, breaks$val)  # no weights
    numKun <- cumsum(wh)
    denKun <- lambda2 * area
    # uncorrected estimate of K
    K <- bind.ratfv(K,
                    data.frame(un=numKun), denKun,
                    "{hat(%s)[%s]^{un}}(r)",
                    "uncorrected estimate of %s",
                    "un",
                    ratio=ratio)
  }
  
  if(any(correction == "border" | correction == "bord.modif")) {
  # border method
  # Compute distances to boundary
    b <- bdist.points(X)
    I <- close$i
    bI <- b[I]
  # apply reduced sample algorithm
    RS <- Kount(DIJ, bI, b, breaks)
    if(any(correction == "bord.modif")) {
      # modified border correction
      denom.area <- eroded.areas(W, r)
      numKbm <- RS$numerator
      denKbm <- lambda2 * denom.area
      K <- bind.ratfv(K,
                      data.frame(bord.modif=numKbm),
                      data.frame(bord.modif=denKbm),
                      "{hat(%s)[%s]^{bordm}}(r)",
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
                      "{hat(%s)[%s]^{bord}}(r)",
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
    denKtrans <- lambda2 * area
    h <- diameter(as.rectangle(W))/2
    numKtrans[r >= h] <- NA
    K <- bind.ratfv(K,
                    data.frame(trans=numKtrans),
                    denKtrans,
                    "{hat(%s)[%s]^{trans}}(r)",
                    "translation-corrected estimate of %s",
                    "trans",
                    ratio=ratio)
  }
  if(any(correction == "isotropic")) {
    ## Ripley isotropic correction
    XI <- ppp(close$xi, close$yi, window=W, check=FALSE)
    edgewt <- edge.Ripley(XI, matrix(DIJ, ncol=1))
    wh <- whist(DIJ, breaks$val, edgewt)
    numKiso <- cumsum(wh)
    denKiso <- lambda2 * area
    h <- diameter(W)/2
    numKiso[r >= h] <- NA
    K <- bind.ratfv(K,
                 data.frame(iso=numKiso),
                 denKiso,
                 "{hat(%s)[%s]^{iso}}(r)",
                 "Ripley isotropic correction estimate of %s",
                 "iso",
                 ratio=ratio)
  }
  #
  # default plot will display all edge corrections
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
