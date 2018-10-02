##
## nnorient.R
##
## nearest neighbour pair orientation distribution
##
## Function \vartheta(phi) defined in
## Illian et al (2008) equ (4.5.3) page 253
##
##  $Revision: 1.4 $ $Date: 2018/10/02 01:21:40 $

nnorient <- function(X, ..., cumulative=FALSE, correction, k = 1,
                     unit=c("degree", "radian"),
                     domain=NULL, ratio=FALSE) {
  stopifnot(is.ppp(X))
  check.1.integer(k)
  stopifnot(k>=1)
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
    correction <- c("bord.modif", "none")

  correction <- pickoption("correction", correction,
                           c(none="none",
                             bord.modif="bord.modif",
                             good="good",
                             best="best"),
                           multi=TRUE)
  correction[correction %in% c("good", "best")] <- "bord.modif"

  ## process point pattern
  Xcoord <- coords(X)
  Ycoord <- Xcoord[nnwhich(X, k=k), ]
  if(!is.null(domain)) {
    inD <- inside.owin(Xcoord$x, Xcoord$y, domain)
    Xcoord <- Xcoord[inD,]
    Ycoord <- Ycoord[inD,]
  } 
  
  dYX <- Ycoord-Xcoord
  ANGLE <- with(dYX, atan2(y, x) * Convert) %% FullCircle
  nangles <- length(ANGLE)
  
  ## initialise output object
  Nphi <- 512
  breaks <- make.even.breaks(bmax=FullCircle, npos=Nphi-1)
  phi <- breaks$r
  Odf <- data.frame(phi  = phi,
                    theo = (if(cumulative) phi else 1)/FullCircle)
  desc <- c("angle argument phi",
            "theoretical isotropic %s")
  NOletter <- if(cumulative) "Theta" else "vartheta"
  NOsymbol <- as.name(NOletter)
  NNO <- ratfv(Odf, NULL, denom=nangles,
              argu="phi",
              ylab=substitute(fn(phi), list(fn=NOsymbol)),
              valu="theo",
              fmla = . ~ phi,
              alim = c(0, FullCircle),
              c("phi",
                "{%s[%s]^{pois}}(phi)"),
              desc,
              fname=NOletter,
              yexp=substitute(fn(phi), list(fn=NOsymbol)))

  ## ^^^^^^^^^^^^^^^  Compute edge corrected estimates ^^^^^^^^^^^^^^^^
  
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
    NNO <- bind.ratfv(NNO,
                     data.frame(un=num.un), den.un,
                    "{hat(%s)[%s]^{un}}(phi)",
                    "uncorrected estimate of %s",
                    "un",
                    ratio=ratio)
  }

  if("bord.modif" %in% correction) {
    ## border type correction
    bX <- bdist.points(X)
    nndX <- nndist(X, k=k)
    if(!is.null(domain)) {
      bX <- bX[inD]
      nndX <- nndX[inD]
    }
    ok <- (nndX < bX)
    nok <- sum(ok)

    rr <- seq(0, max(bX), length=256)

    if(nok == 0) {
      num.bm <- numeric(Nphi) # i.e. rep(0, Nphi)
    } else {
      Ar <- eroded.areas(W, rr)
      Arf <- approxfun(rr, Ar, rule=2)
      AI <- Arf(bX)
      edgewt <- ifelse(ok, pmin(area(W)/AI, 100), 0)
      if(cumulative) {
        wh <- whist(ANGLE, breaks$val, edgewt)
        num.bm <- cumsum(wh)/mean(edgewt)
      } else {
        w <- edgewt/sum(edgewt)
        kd <- circdensity(ANGLE, ..., weights=w, n=Nphi, unit=unit)
        num.bm <- kd$y * nok
      }
    }
    den.bm <- nok
    NNO <- bind.ratfv(NNO,
                      data.frame(bordm=num.bm),
                      den.bm,
                      "{hat(%s)[%s]^{bordm}}(phi)",
                      "modified border-corrected estimate of %s",
                      "bordm",
                      ratio=ratio)
  }
 
  unitname(NNO) <- switch(unit,
                         degree = c("degree", "degrees"),
                         radian = c("radian", "radians"))
  return(NNO)
}
