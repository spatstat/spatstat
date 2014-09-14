##
## ppod.R
##
## point pair orientation distribution
##
## Function O_{r1,r2}(phi) defined in
## Stoyan & Stoyan (1994) equ (14.53) page 271
##

pairorient <- function(X, r1, r2, ...,
                       correction, ratio=FALSE,
                       units=c("degrees", "radians"),
                       domain=NULL) {
  stopifnot(is.ppp(X))
  check.1.real(r1)
  check.1.real(r2)
  stopifnot(r1 < r2)
  W <- Window(X)
  if(!is.null(domain))
    stopifnot(is.subset.owin(domain, Window(X)))
  
  units <- match.arg(units)
  switch(units,
         degrees = {
           FullCircle <- 360
           Convert <- 180/pi
         },
         radians = {
           FullCircle <- 2 * pi
           Convert <- 1
         })

  ## choose correction(s)
  correction.given <- !missing(correction) && !is.null(correction)
  if(!correction.given)
    correction <- c("isotropic", "translate")
  correction <- pickoption("correction", correction,
                           c(none="none",
##                             border="border",
##                             "bord.modif"="bord.modif",
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

  
  ## border corrections not implemented
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
  breaks <- make.even.breaks(bmax=FullCircle, npos=511)
  phi <- breaks$r
  Odf <- data.frame(phi  = phi,
                    theo = phi/FullCircle)
  desc <- c("distance argument r",
            "theoretical isotropic %s")
  OO <- ratfv(Odf, NULL, denom=nrow(close),
              argu="phi",
              ylab=substitute(O[R1,R2](r), list(R1=r1, R2=r2)), 
              valu="theo",
              fmla = . ~ phi,
              alim = c(0, FullCircle),
              c("phi",
                "{%s[%s]^{pois}}(phi)"),
              desc,
              fname=c("O", paste0("list(", r1, ",", r2, ")")),
              yexp=substitute(O[list(R1,R2)](phi),
                list(R1=r1,R2=r2)))

  ## ^^^^^^^^^^^^^^^  Compute edge corrected estimates ^^^^^^^^^^^^^^^^
  
  if(any(correction == "none")) {
    ## uncorrected! For demonstration purposes only!
    wh <- whist(ANGLE, breaks$val)  # no weights
    num.un <- cumsum(wh)
    den.un <- length(ANGLE)
    ## uncorrected estimate 
    OO <- bind.ratfv(OO,
                     data.frame(un=num.un), den.un,
                    "{hat(%s)[%s]^{un}}(phi)",
                    "uncorrected estimate of %s",
                    "un",
                    ratio=ratio)
  }

  if(any(correction == "translate")) {
    ## Ohser-Stoyan translation correction
    edgewt <- edge.Trans(dx=close$dx, dy=close$dy, W=W, paired=TRUE)
    wh <- whist(ANGLE, breaks$val, edgewt)
    num.trans <- cumsum(wh)/mean(edgewt)
    den.trans <- length(ANGLE)
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
    wh <- whist(ANGLE, breaks$val, edgewt)
    num.iso <- cumsum(wh)/mean(edgewt)
    den.iso <- length(ANGLE)
    OO <- bind.ratfv(OO,
                     data.frame(iso=num.iso),
                     den.iso,
                     "{hat(%s)[%s]^{iso}}(phi)",
                     "Ripley isotropic-corrected estimate of %s",
                     "iso",
                     ratio=ratio)
  }
  unitname(OO) <- switch(units,
                         degrees = c("degree", "degrees"),
                         radians = c("radian", "radians"))
  return(OO)
}


