#
# linearK
#
# $Revision: 1.30 $ $Date: 2013/04/25 06:37:43 $
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
  sX <- summary(X)
  np <- sX$npoints
  lengthL <- sX$totlength
  # compute K
  denom <- np * (np - 1)/lengthL
  K <- linearKengine(X, r=r, denom=denom, correction=correction, ...)
  # set appropriate y axis label
  switch(correction,
         Ang  = {
           ylab <- quote(K[L](r))
           fname <- "K[L]"
         },
         none = {
           ylab <- quote(K[net](r))
           fname <- "K[net]"
         })
  K <- rebadge.fv(K, new.ylab=ylab, new.fname=fname)
  return(K)
}

linearKinhom <- function(X, lambda=NULL, r=NULL,  ...,
                         correction="Ang", normalise=TRUE) {
  stopifnot(inherits(X, "lpp"))
  correction <- pickoption("correction", correction,
                           c(none="none",
                             Ang="Ang",
                             best="Ang"),
                           multi=FALSE)
  if(is.null(lambda))
    linearK(X, r=r, ..., correction=correction)
  # extract info about pattern
  sX <- summary(X)
  np <- sX$npoints
  lengthL <- sX$totlength
  #
  lambdaX <- getlambda.lpp(lambda, X, ...)
  #
  invlam <- 1/lambdaX
  invlam2 <- outer(invlam, invlam, "*")
  denom <- if(!normalise) lengthL else sum(invlam)
  K <- linearKengine(X, ...,
                     r=r, reweight=invlam2, denom=denom, correction=correction)
  # set appropriate y axis label
  switch(correction,
         Ang  = {
           ylab <- quote(K[LI](r))
           fname <- "K[LI]"
         },
         none = {
           ylab <- quote(K[netI](r))
           fname <- "K[netI]"
         })
  K <- rebadge.fv(K, new.fname=fname, new.ylab=ylab)
  return(K)
}

getlambda.lpp <- function(lambda, X, ...) {
  lambdaname <- deparse(substitute(lambda))
  XX <- as.ppp(X)
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

  return(lambdaX)
}

linearKengine <- function(X, ..., r=NULL, reweight=NULL, denom=1,
                          correction="Ang", showworking=FALSE) {
  # extract info about pattern
  sX <- summary(X)
  np <- sX$npoints
  lengthL <- sX$totlength
  # extract linear network
  L <- X$domain
  # extract points
  Y <- as.ppp(X)
  W <- Y$window
  # determine r values
  rmaxdefault <- 0.98 * circumradius(L)
  breaks <- handle.r.b.args(r, NULL, W, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max
  #
  if(np < 2) {
    # no pairs to count: return zero function
    zeroes <- numeric(length(r))
    df <- data.frame(r = r, est = zeroes)
    K <- fv(df, "r", substitute(linearK(r), NULL),
            "est", . ~ r, c(0, rmax),
            c("r", "%s(r)"),
            c("distance argument r", "estimated %s"),
            fname = "linearK")
    return(K)
  }
  # compute pairwise distances  
  D <- pairdist(X)
  #---  compile into K function ---
  if(correction == "none" && is.null(reweight)) {
    # no weights (Okabe-Yamada)
    K <- compileK(D, r, denom=denom)
    unitname(K) <- unitname(X)
    return(K)
  }
  if(correction == "none")
     edgewt <- 1
  else {
     # inverse m weights (Wei's correction)
     # compute m[i,j]
     m <- matrix(1, np, np)
     for(j in 1:np) 
       m[ -j, j] <- countends(L, Y[-j], D[-j,j])
     if(any(uhoh <- (m == 0))) {
       warning("Internal error: disc boundary count equal to zero")
       m[uhoh] <- 1
     }
     edgewt <- 1/m
  }
  # compute K
  wt <- if(!is.null(reweight)) edgewt * reweight else edgewt
  K <- compileK(D, r, weights=wt, denom=denom)
  # tack on theoretical value
  K <- bind.fv(K, data.frame(theo=r), "%s[theo](r)",
               "theoretical Poisson %s")
  unitname(K) <- unitname(X)
  fvnames(K, ".") <- rev(fvnames(K, "."))
  # show working
  if(showworking)
    attr(K, "working") <- list(D=D, wt=wt)
  return(K)
}

