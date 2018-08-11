##
## pairorient.R
##
## point pair orientation distribution
##
## Function O_{r1,r2}(phi) defined in
## Stoyan & Stoyan (1994) equ (14.53) page 271
##
##     and its derivative estimated by kernel smoothing
##
##  $Revision: 1.9 $ $Date: 2014/12/05 06:59:53 $

pairorient <- function(X, r1, r2, ...,
                       cumulative=FALSE,
                       correction, ratio=FALSE,
                       unit=c("degree", "radian"),
                       domain=NULL) {
  stopifnot(is.ppp(X))
  check.1.real(r1)
  check.1.real(r2)
  stopifnot(r1 < r2)
  W <- Window(X)
  if(!is.null(domain))
    stopifnot(is.subset.owin(domain, W))
  
  unit <- match.arg(unit)
  switch(unit,
         degree = {
           FullCircle <- 360
           Convert <- 180/pi
         },
         radian = {
           FullCircle <- 2 * pi
           Convert <- 1
         })

  ## choose correction(s)
  correction.given <- !missing(correction) && !is.null(correction)
  if(!correction.given)
    correction <- c("border", "isotropic", "translate")
  correction <- pickoption("correction", correction,
                           c(none="none",
                             border="border",
                             bord.modif="bord.modif",
                             isotropic="isotropic",
                             Ripley="isotropic",
                             trans="translate",
                             translate="translate",
                             translation="translate",
                             good="good",
                             best="best"),
                           multi=TRUE)
#  best.wanted <- ("best" %in% correction)
  ## replace 'good' by the optimal choice for this size of dataset
  if("good" %in% correction)
    correction[correction == "good"] <- good.correction.K(X)
  ## retain only corrections that are implemented for the window
  correction <- implemented.for.K(correction, W$type, correction.given)

  
  ## Find close pairs in range [r1, r2]
  close <- as.data.frame(closepairs(X, r2))
  ok <- with(close, r1 <= d & d <= r2)
  if(!is.null(domain))
      ok <- ok & with(close, inside.owin(xi, yi, domain))
  if(!any(ok)) {
    warning(paste("There are no pairs of points in the distance range",
                  prange(c(r1,r2))))
    return(NULL)
  }
  close <- close[ok, , drop=FALSE]
  ANGLE <- with(close, atan2(dy, dx) * Convert) %% FullCircle

  ## initialise output object
  Nphi <- 512
  breaks <- make.even.breaks(bmax=FullCircle, npos=Nphi-1)
  phi <- breaks$r
  Odf <- data.frame(phi  = phi,
                    theo = (if(cumulative) phi else 1)/FullCircle)
  desc <- c("angle argument phi",
            "theoretical isotropic %s")
  Oletter <- if(cumulative) "O" else "o"
  Osymbol <- as.name(Oletter)
  OO <- ratfv(Odf, NULL, denom=nrow(close),
              argu="phi",
              ylab=substitute(fn[R1,R2](phi), list(R1=r1, R2=r2, fn=Osymbol)),
              valu="theo",
              fmla = . ~ phi,
              alim = c(0, FullCircle),
              c("phi",
                "{%s[%s]^{pois}}(phi)"),
              desc,
              fname=c(Oletter, paste0("list(", r1, ",", r2, ")")),
              yexp=substitute(fn[list(R1,R2)](phi),
                list(R1=r1,R2=r2,fn=Osymbol)))

  ## ^^^^^^^^^^^^^^^  Compute edge corrected estimates ^^^^^^^^^^^^^^^^

  nangles <- length(ANGLE)
  
  if(any(correction == "none")) {
    ## uncorrected! For demonstration purposes only!
    if(cumulative) {
      wh <- whist(ANGLE, breaks$val)  # no weights
      num.un <- cumsum(wh)
    } else {
      kd <- circdensity(ANGLE, ..., n=Nphi, unit=unit)
      num.un <- kd$y * nangles
    }
    den.un <- nangles
    ## uncorrected estimate 
    OO <- bind.ratfv(OO,
                     data.frame(un=num.un), den.un,
                    "{hat(%s)[%s]^{un}}(phi)",
                    "uncorrected estimate of %s",
                    "un",
                    ratio=ratio)
  }

  if(any(c("border", "bord.modif") %in% correction)) {
    ## border type corrections
    bX <- bdist.points(X)
    bI <- bX[close$i]
    if("border" %in% correction) {
      bok <- (bI > r2)
      ANGLEok <- ANGLE[bok]
      nok <- length(ANGLEok)
      if(cumulative) {
        wh <- whist(ANGLEok, breaks$val)
        num.bord <- cumsum(wh)
      } else {
        kd <- circdensity(ANGLEok, ..., n=Nphi, unit=unit)
        num.bord <- kd$y * nok
      }
      den.bord <- nok
      OO <- bind.ratfv(OO,
                       data.frame(border=num.bord),
                       den.bord,
                       "{hat(%s)[%s]^{bord}}(phi)",
                       "border-corrected estimate of %s",
                       "border",
                       ratio=ratio)
    }
    if("bord.modif" %in% correction) {
      ok <- (close$d < bI)
      nok <- sum(ok)
      inradius <- max(distmap(W, invert=TRUE))
      rrr <- range(r2, inradius)
      rr <- seq(rrr[1], rrr[2], length=256)
      Ar <- eroded.areas(W, rr)
      Arf <- approxfun(rr, Ar, rule=2)
      AI <- (Arf(bX))[close$i]
      edgewt <- ifelse(ok, pmin(area(W)/AI, 100), 0)
      if(cumulative) {
        wh <- whist(ANGLE, breaks$val, edgewt)
        num.bm <- cumsum(wh)/mean(edgewt)
      } else {
        w <- edgewt/sum(edgewt)
        kd <- circdensity(ANGLE, ..., weights=w, n=Nphi, unit=unit)
        num.bm <- kd$y * nok
      }
      den.bm <- nok
      OO <- bind.ratfv(OO,
                       data.frame(bordm=num.bm),
                       den.bm,
                       "{hat(%s)[%s]^{bordm}}(phi)",
                       "modified border-corrected estimate of %s",
                       "bordm",
                       ratio=ratio)
    }
  }
  if(any(correction == "translate")) {
    ## Ohser-Stoyan translation correction
    edgewt <- edge.Trans(dx=close$dx, dy=close$dy, W=W, paired=TRUE)
    if(cumulative) {
      wh <- whist(ANGLE, breaks$val, edgewt)
      num.trans <- cumsum(wh)/mean(edgewt)
    } else {
      w <- edgewt/sum(edgewt)
      kd <- circdensity(ANGLE, ..., weights=w, n=Nphi, unit=unit)
      num.trans <- kd$y * nangles
    }
    den.trans <- nangles
    OO <- bind.ratfv(OO,
                     data.frame(trans=num.trans),
                     den.trans,
                     "{hat(%s)[%s]^{trans}}(phi)",
                     "translation-corrected estimate of %s",
                     "trans",
                     ratio=ratio)
  }
  if(any(correction == "isotropic")) {
    ## Ripley isotropic correction
    XI <- ppp(close$xi, close$yi, window=W, check=FALSE)
    DIJ <- close$d
    edgewt <- edge.Ripley(XI, matrix(DIJ, ncol=1))
    if(cumulative) {
      wh <- whist(ANGLE, breaks$val, edgewt)
      num.iso <- cumsum(wh)/mean(edgewt)
    } else {
      w <- edgewt/sum(edgewt)
      kd <- circdensity(ANGLE, ..., weights=w, n=Nphi, unit=unit)
      num.iso <- kd$y * nangles
    }
    den.iso <- nangles
    OO <- bind.ratfv(OO,
                     data.frame(iso=num.iso),
                     den.iso,
                     "{hat(%s)[%s]^{iso}}(phi)",
                     "Ripley isotropic-corrected estimate of %s",
                     "iso",
                     ratio=ratio)
  }
  unitname(OO) <- switch(unit,
                         degree = c("degree", "degrees"),
                         radian = c("radian", "radians"))
  return(OO)
}
