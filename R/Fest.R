#
#	Fest.S
#
#	S function empty.space()
#	Computes estimates of the empty space function
#
#	$Revision: 4.32 $	$Date: 2013/04/25 06:37:43 $
#
"Fest" <- 	
"empty.space" <-
function(X, ..., eps = NULL, r=NULL, breaks=NULL,
         correction=c("rs", "km", "cs")) {
  verifyclass(X, "ppp")
  
# Intensity estimate
  W <- X$window
  npts <- npoints(X)
  lambda <- npts/area.owin(W)
  
# First discretise
  dwin <- as.mask(W, eps)
  dX <- ppp(X$x, X$y, window=dwin, check=FALSE)
#        
# histogram breakpoints 
#
  rmaxdefault <- rmax.rule("F", dwin, lambda)
  breaks <- handle.r.b.args(r, breaks, dwin, eps,
                                  rmaxdefault=rmaxdefault)
  rvals <- breaks$r
  rmax  <- breaks$max
  
# choose correction(s)
  correction.given <- !missing(correction) && !is.null(correction)
  if(is.null(correction)) {
    correction <- c("rs", "km", "cs")
  } else correction <- pickoption("correction", correction,
                           c(none="none",
                             border="rs",
                             rs="rs",
                             KM="km",
                             km="km",
                             Kaplan="km",
                             cs="cs",
                             ChiuStoyan="cs",
                             Hanisch="cs",
                             han="cs",
                             best="km"),
                           multi=TRUE)
  
# initialise fv object
  df <- data.frame(r=rvals, theo=1-exp(-lambda * pi * rvals^2))
  Z <- fv(df, "r", substitute(F(r), NULL), "theo", . ~ r,
          c(0,rmax),
          c("r", "%s[pois](r)"), 
          c("distance argument r", "theoretical Poisson %s"),
          fname="F")
  nr <- length(rvals)
  zeroes <- numeric(nr)
#
#  compute distances and censoring distances
  if(X$window$type == "rectangle") {
    # original data were in a rectangle
    # output of exactdt() is sufficient
    e <- exactdt(dX)
    dist <- e$d
    bdry <- e$b
  } else {
    # window is irregular..
    # Distance transform & boundary distance for all pixels
    e <- exactdt(dX)
    b <- bdist.pixels(dX$window, style="matrix")
    # select only those pixels inside mask
    mm <- dwin$m
    dist <- e$d[mm]
    bdry <- b[mm]
  }
  
# censoring indicators
  d <- (dist <= bdry)
#  observed distances
  o <- pmin(dist, bdry)

### start calculating estimates of F
  
  if("none" %in% correction) {
    #  UNCORRECTED e.d.f. of empty space distances
    if(npts == 0)
      edf <- zeroes
    else {
      hh <- hist(dist[dist <= rmax],breaks=breaks$val,plot=FALSE)$counts
      edf <- cumsum(hh)/length(dist)
    }
    Z <- bind.fv(Z, data.frame(raw=edf), "hat(%s)[raw](r)",
                 "uncorrected estimate of %s", "raw")
  }
  
  if("cs" %in% correction) {
    # Chiu-Stoyan correction
    if(npts == 0)
      cs <- zeroes
    else {
      #  uncensored distances
      x <- dist[d]
      #  weights
      a <- eroded.areas(W, rvals)
      # calculate Hanisch estimator
      h <- hist(x[x <= rmax], breaks=breaks$val, plot=FALSE)$counts
      H <- cumsum(h/a)
      cs <- H/max(H[is.finite(H)])
    }
    # add to fv object
    Z <- bind.fv(Z, data.frame(cs=cs),
                 "hat(%s)[cs](r)", 
                 "Chiu-Stoyan estimate of %s",
                 "cs")
  }

  if(any(correction %in% c("rs", "km"))) {
    # calculate Kaplan-Meier and/or border corrected (Reduced Sample) estimators
    want.rs <- "rs" %in% correction
    want.km <- "km" %in% correction
    selection <- c(want.rs, want.km, want.km)
    tags <- c("rs", "km", "hazard")[selection]
    labels <- c("hat(%s)[bord](r)", "hat(%s)[km](r)", "hazard(r)")[selection]
    descr <- c("border corrected estimate of %s",
               "Kaplan-Meier estimate of %s",
               "Kaplan-Meier estimate of hazard function lambda(r)")[selection]
    if(npts == 0) {
      result <- as.data.frame(matrix(0, nr, length(tags)))
      names(result) <- tags
    } else {
      result <- km.rs.opt(o, bdry, d, breaks, KM=want.km, RS=want.rs)
      result <- as.data.frame(result[tags])
    }
    # add to fv object
    Z <- bind.fv(Z, result,
                 labels, descr, if(want.km) "km" else "rs")
  }
  
  # wrap up
  unitname(Z) <- unitname(X)
  
  # remove 'hazard' from the dotnames
  nama <- names(Z)
  fvnames(Z, ".") <- rev(nama[!(nama %in% c("r", "hazard"))])
  
  # determine recommended plot range
  attr(Z, "alim") <- with(Z, range(.x[is.finite(.y) & .y <= 0.9]))

  return(Z)
}

	
