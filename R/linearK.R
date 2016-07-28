#
# linearK
#
# $Revision: 1.40 $ $Date: 2016/07/18 03:48:10 $
#
# K function for point pattern on linear network
#
#
linearK <- function(X, r=NULL, ..., correction="Ang") {
  stopifnot(inherits(X, "lpp"))
  correction <- pickoption("correction", correction,
                           c(none="none",
                             Ang="Ang",
                             best="Ang"),
                           multi=FALSE)
  # extract info about pattern
  np <- npoints(X)
  lengthL <- volume(domain(X))
  # compute K
  denom <- np * (np - 1)/lengthL
  K <- linearKengine(X, r=r, denom=denom, correction=correction, ...)
  # set appropriate y axis label
  switch(correction,
         Ang  = {
           ylab <- quote(K[L](r))
           fname <- c("K", "L")
         },
         none = {
           ylab <- quote(K[net](r))
           fname <- c("K", "net")
         })
  K <- rebadge.fv(K, new.ylab=ylab, new.fname=fname)
  return(K)
}

linearKinhom <- function(X, lambda=NULL, r=NULL,  ...,
                         correction="Ang", normalise=TRUE, normpower=1) {
  stopifnot(inherits(X, "lpp"))
  correction <- pickoption("correction", correction,
                           c(none="none",
                             Ang="Ang",
                             best="Ang"),
                           multi=FALSE)
  if(is.null(lambda))
    linearK(X, r=r, ..., correction=correction)
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
  K <- linearKengine(X, ...,
                     r=r, reweight=invlam2, denom=denom, correction=correction)
  # set appropriate y axis label
  switch(correction,
         Ang  = {
           ylab <- quote(K[L, inhom](r))
           fname <- c("K", "list(L, inhom)")
         },
         none = {
           ylab <- quote(K[net, inhom](r))
           fname <- c("K", "list(net, inhom)")
         })
  K <- rebadge.fv(K, new.fname=fname, new.ylab=ylab)
  attr(K, "dangerous") <- attr(lambdaX, "dangerous")
  return(K)
}

getlambda.lpp <- function(lambda, X, ..., update=TRUE) {
  lambdaname <- deparse(substitute(lambda))
  XX <- as.ppp(X)
  missup <- missing(update)
  if(update && (is.lppm(lambda) || is.ppm(lambda))) {
    danger <- FALSE
    lambda <- if(is.lppm(lambda)) update(lambda, X) else update(lambda, XX)
    if(missup)
      warn.once("lin.inhom.update",
                "The behaviour of linearKinhom and similar functions",
                "when lambda is an lppm object",
                "has changed (in spatstat 1.41-0 and later)",
                "See help(linearKinhom)")
  } else danger <- TRUE
  lambdaX <-
    if(is.vector(lambda)) lambda  else
    if(is.function(lambda)) lambda(XX$x, XX$y, ...) else
    if(is.im(lambda)) safelookup(lambda, XX) else 
    if(inherits(lambda, "linim")) safelookup(as.im(lambda), XX) else 
    if(is.ppm(lambda) || inherits(lambda, "lppm"))
      predict(lambda, locations=as.data.frame(XX)) else
    stop(paste(lambdaname, "should be",
               "a numeric vector, function, pixel image, or fitted model"))
  if(!is.numeric(lambdaX))
    stop(paste("Values of", lambdaname, "are not numeric"))
  if((nv <- length(lambdaX)) != (np <- npoints(X)))
     stop(paste("Obtained", nv, "values of", lambdaname,
	   "but point pattern contains", np, "points"))
  if(any(lambdaX < 0))
    stop(paste("Negative values of", lambdaname, "obtained"))
  if(any(lambdaX == 0))
    stop(paste("Zero values of", lambdaname, "obtained"))
  if(danger)
    attr(lambdaX, "dangerous") <- lambdaname
  return(lambdaX)
}

linearKengine <- function(X, ..., r=NULL, reweight=NULL, denom=1,
                          correction="Ang", showworking=FALSE) {
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
  fname <- c("K", type)
  ylab <- substitute(K[type](r), list(type=type))
  #
  if(np < 2) {
    # no pairs to count: return zero function
    zeroes <- numeric(length(r))
    df <- data.frame(r = r, est = zeroes)
    K <- fv(df, "r", ylab,
            "est", . ~ r, c(0, rmax),
            c("r", makefvlabel(NULL, "hat", fname)), 
            c("distance argument r", "estimated %s"),
            fname = fname)
    return(K)
  }
  # compute pairwise distances  
  D <- pairdist(X)
  #---  compile into K function ---
  if(correction == "none" && is.null(reweight)) {
    # no weights (Okabe-Yamada)
    K <- compileK(D, r, denom=denom, fname=fname)
    K <- rebadge.fv(K, ylab, fname)
    unitname(K) <- unitname(X)
    return(K)
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
     if(any(uhoh <- (m == 0))) {
       warning("Internal error: disc boundary count equal to zero")
       m[uhoh] <- 1
     }
     edgewt <- 1/m
  }
  # compute K
  wt <- if(!is.null(reweight)) edgewt * reweight else edgewt
  K <- compileK(D, r, weights=wt, denom=denom, fname=fname)
  # tack on theoretical value
  K <- bind.fv(K, data.frame(theo=r),
               makefvlabel(NULL, NULL, fname, "theo"),
               "theoretical Poisson %s")
  K <- rebadge.fv(K, ylab, fname)
  unitname(K) <- unitname(X)
  fvnames(K, ".") <- rev(fvnames(K, "."))
  # show working
  if(showworking)
    attr(K, "working") <- list(D=D, wt=wt)
  attr(K, "correction") <- correction
  return(K)
}

