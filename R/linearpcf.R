#
# linearpcf.R
#
# $Revision: 1.18 $ $Date: 2016/07/18 03:48:05 $
#
# pair correlation function for point pattern on linear network
#
#
linearpcf <- function(X, r=NULL, ..., correction="Ang") {
  stopifnot(inherits(X, "lpp"))
  correction <- pickoption("correction", correction,
                           c(none="none",
                             Ang="Ang",
                             best="Ang"),
                           multi=FALSE)
  # extract info about pattern
  np <- npoints(X)
  lengthL <- volume(domain(X))
  # compute
  denom <- np * (np - 1)/lengthL
  g <- linearpcfengine(X, r=r, ..., denom=denom, correction=correction)
  # set appropriate y axis label
  switch(correction,
         Ang  = {
           ylab <- quote(g[L](r))
           fname <- c("g", "L")
         },
         none = {
           ylab <- quote(g[net](r))
           fname <- c("g", "net")
         })
  g <- rebadge.fv(g, new.ylab=ylab, new.fname=fname)
  return(g)
}

linearpcfinhom <- function(X, lambda=NULL, r=NULL,  ...,
                           correction="Ang", normalise=TRUE, normpower=1) {
  stopifnot(inherits(X, "lpp"))
  correction <- pickoption("correction", correction,
                           c(none="none",
                             Ang="Ang",
                             best="Ang"),
                           multi=FALSE)
  if(is.null(lambda))
    linearpcf(X, r=r, ..., correction=correction)
  if(normalise) {
    check.1.real(normpower)
    stopifnot(normpower >= 1)
  }
  # extract info about pattern
  lengthL <- volume(domain(X))
  #
  lambdaX <- getlambda.lpp(lambda, X, ...)
  #
  invlam <- 1/lambdaX
  invlam2 <- outer(invlam, invlam, "*")
  denom <- if(!normalise) lengthL else
           if(normpower == 1) sum(invlam) else
           lengthL * (sum(invlam)/lengthL)^normpower
  g <- linearpcfengine(X, ..., r=r,
                       reweight=invlam2, denom=denom, correction=correction)
  # extract bandwidth
  bw <- attr(g, "bw")
  # set appropriate y axis label
  switch(correction,
         Ang  = {
           ylab <- quote(g[L, inhom](r))
           fname <- c("g", "list(L, inhom)")
         },
         none = {
           ylab <- quote(g[net, inhom](r))
           fname <- c("g", "list(net, inhom)")
         })
  g <- rebadge.fv(g, new.fname=fname, new.ylab=ylab)
  # reattach bandwidth
  attr(g, "bw") <- bw
  attr(g, "dangerous") <- attr(lambdaX, "dangerous")
  return(g)
}


linearpcfengine <- function(X, ..., r=NULL,
                            reweight=NULL, denom=1, correction="Ang") {
  # ensure distance information is present
  X <- as.lpp(X, sparse=FALSE)
  # extract info about pattern
  np <- npoints(X)
  # extract linear network
  L <- domain(X)
  W <- Window(L)
  # determine r values
  rmaxdefault <- 0.98 * boundingradius(L)
  breaks <- handle.r.b.args(r, NULL, W, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max
  #
  type <- if(correction == "Ang") "L" else "net"
  fname <- c("g", type)
  ylab <- substitute(g[type](r), list(type=type))
  #  
  if(np < 2) {
    # no pairs to count: return zero function
    zeroes <- numeric(length(r))
    df <- data.frame(r = r, est = zeroes)
    g <- fv(df, "r", ylab,
            "est", . ~ r, c(0, rmax),
            c("r", makefvlabel(NULL, "hat", fname)), 
            c("distance argument r", "estimated %s"),
            fname = fname)
    return(g)
  }
  # compute pairwise distances  
  D <- pairdist(X)
  #---  compile into pcf ---
  if(correction == "none" && is.null(reweight)) {
    # no weights (Okabe-Yamada)
    g <- compilepcf(D, r, denom=denom, fname=fname)
    unitname(g) <- unitname(X)
    attr(g, "correction") <- correction
    return(g)
  }
  if(correction == "none")
     edgewt <- 1
  else {
     # inverse m weights (Wei's correction)
     # determine tolerance
     toler <- default.linnet.tolerance(L)
     # compute m[i,j]
     m <- matrix(1, np, np)
     for(j in 1:np) 
       m[ -j, j] <- countends(L, X[-j], D[-j,j], toler=toler)
     edgewt <- 1/m
  }
  # compute pcf
  wt <- if(!is.null(reweight)) edgewt * reweight else edgewt
  g <- compilepcf(D, r, weights=wt, denom=denom, ..., fname=fname)
  # extract bandwidth
  bw <- attr(g, "bw")
  # tack on theoretical value
  g <- bind.fv(g, data.frame(theo=rep.int(1,length(r))),
               makefvlabel(NULL, NULL, fname, "pois"),
               "theoretical Poisson %s")
  # tweak
  unitname(g) <- unitname(X)
  fvnames(g, ".") <- rev(fvnames(g, "."))
  # tack on bandwidth again
  attr(g, "bw") <- bw
  attr(g, "correction") <- correction
  return(g)
}

