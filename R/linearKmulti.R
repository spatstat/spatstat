#
# linearKmulti
#
# $Revision: 1.2 $ $Date: 2014/02/16 08:50:46 $
#
# K functions for multitype point pattern on linear network
#
#

linearKdot <- function(X, i, r=NULL, ..., correction="Ang") {
  if(!is.multitype(X, dfok=FALSE)) 
	stop("Point pattern must be multitype")
  marx <- marks(X)
  lev <- levels(marx)
  if(missing(i)) i <- lev[1] else
    if(!(i %in% lev)) stop(paste("i = ", i , "is not a valid mark"))  
  I <- (marx == i)
  J <- rep(TRUE, npoints(X))  # i.e. all points
  result <- linearKmulti(X, I, J,
                         r=r, correction=correction, ...)
  correction <- attr(result, "correction")
  type <- if(correction == "Ang") "L" else "net"
  result <- rebadge.as.dotfun(result, "K", type, i)
  return(result)
}

linearKcross <- function(X, i, j, r=NULL, ..., correction="Ang") {
  if(!is.multitype(X, dfok=FALSE)) 
	stop("Point pattern must be multitype")
  marx <- marks(X)
  lev <- levels(marx)
  if(missing(i)) i <- lev[1] else
    if(!(i %in% lev)) stop(paste("i = ", i , "is not a valid mark"))
  if(missing(j)) j <- lev[2] else
    if(!(j %in% lev)) stop(paste("j = ", j , "is not a valid mark"))
  #
  if(i == j) {
    result <- linearK(X[marx == i], r=r, correction=correction, ...)
  } else {
    I <- (marx == i)
    J <- (marx == j)
    result <- linearKmulti(X, I, J, r=r, correction=correction, ...)
  }
  # rebrand
  correction <- attr(result, "correction")
  type <- if(correction == "Ang") "L" else "net"
  result <- rebadge.as.crossfun(result, "K", type, i, j)
  return(result)
}

linearKmulti <- function(X, I, J, r=NULL, ..., correction="Ang") {
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
  # validate I, J
  if(!is.logical(I) || !is.logical(J))
    stop("I and J must be logical vectors")
  if(length(I) != np || length(J) != np)
    stop(paste("The length of I and J must equal",
               "the number of points in the pattern"))
	
  if(!any(I)) stop("no points satisfy I")
#  if(!any(J)) stop("no points satisfy J")
		
  nI <- sum(I)
  nJ <- sum(J)
  nIandJ <- sum(I & J)
  lambdaI <- nI/lengthL
  lambdaJ <- nJ/lengthL
  # compute K
  denom <- (nI * nJ - nIandJ)/lengthL
  K <- linearKmultiEngine(X, I, J, r=r, denom=denom,
                          correction=correction, ...)
  # set appropriate y axis label
  correction <- attr(K, "correction")
  type <- if(correction == "Ang") "L" else "net"
  K <- rebadge.as.crossfun(K, "K", type, "I", "J")
  return(K)
}

# ................ inhomogeneous ............................

linearKdot.inhom <- function(X, i, lambdaI, lambdadot,
                             r=NULL, ..., correction="Ang", normalise=TRUE) {
  if(!is.multitype(X, dfok=FALSE)) 
	stop("Point pattern must be multitype")
  marx <- marks(X)
  lev <- levels(marx)
  if(missing(i)) i <- lev[1] else
    if(!(i %in% lev)) stop(paste("i = ", i , "is not a valid mark"))  
  I <- (marx == i)
  J <- rep(TRUE, npoints(X))  # i.e. all points
  # for better error messages
  lambdadot <- getlambda.lpp(lambdadot, X, ...)
  # compute
  result <- linearKmulti.inhom(X, I, J, lambdaI, lambdadot, 
                               r=r, correction=correction, normalise=normalise,
                               ...)
  ## relabel
  correction <- attr(result, "correction")
  type <- if(correction == "Ang") "L, inhom" else "net, inhom"
  result <- rebadge.as.dotfun(result, "K", type, i)
  return(result)
}

linearKcross.inhom <- function(X, i, j, lambdaI, lambdaJ,
                               r=NULL, ...,
                               correction="Ang", normalise=TRUE) {
  if(!is.multitype(X, dfok=FALSE)) 
	stop("Point pattern must be multitype")
  marx <- marks(X)
  lev <- levels(marx)
  if(missing(i)) i <- lev[1] else
    if(!(i %in% lev)) stop(paste("i = ", i , "is not a valid mark"))
  if(missing(j)) j <- lev[2] else
    if(!(j %in% lev)) stop(paste("j = ", j , "is not a valid mark"))
  #
  if(i == j) {
    I <- (marx == i)
    result <- linearKinhom(X[I], lambda=lambdaI, r=r,
                           correction=correction, normalise=normalise, ...)
  } else {
    I <- (marx == i)
    J <- (marx == j)
    result <- linearKmulti.inhom(X, I, J, lambdaI, lambdaJ,
                                 r=r, correction=correction,
                                 normalise=normalise, ...)
  }
  # rebrand
  correction <- attr(result, "correction")
  type <- if(correction == "Ang") "L, inhom" else "net, inhom"
  result <- rebadge.as.crossfun(result, "K", type, i, j)
  return(result)
}

linearKmulti.inhom <- function(X, I, J, lambdaI, lambdaJ,
                               r=NULL, ...,
                               correction="Ang", normalise=TRUE) {
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
  # validate I, J
  if(!is.logical(I) || !is.logical(J))
    stop("I and J must be logical vectors")
  if(length(I) != np || length(J) != np)
    stop(paste("The length of I and J must equal",
               "the number of points in the pattern"))
	
  if(!any(I)) stop("no points satisfy I")

  # validate lambda vectors
  lambdaI <- getlambda.lpp(lambdaI, X[I], ...)
  lambdaJ <- getlambda.lpp(lambdaJ, X[J], ...)

  # compute K
  weightsIJ <- outer(1/lambdaI, 1/lambdaJ, "*")
  denom <- if(!normalise) lengthL else sum(1/lambdaI)
  K <- linearKmultiEngine(X, I, J, r=r,
                          reweight=weightsIJ, denom=denom,
                          correction=correction, ...)
  # set appropriate y axis label
  correction <- attr(K, "correction")
  type <- if(correction == "Ang") "L, inhom" else "net, inhom"
  K <- rebadge.as.crossfun(K, "K", type, "I", "J")
  return(K)
}

# .............. internal ...............................

linearKmultiEngine <- function(X, I, J, ..., r=NULL, reweight=NULL, denom=1,
                          correction="Ang", showworking=FALSE) {
  # extract info about pattern
  sX <- summary(X)
  np <- sX$npoints
  lengthL <- sX$totlength
  # extract linear network
  L <- X$domain
  # extract points
  XP <- as.ppp(X)
  W <- as.owin(XP)
  # determine r values
  rmaxdefault <- 0.98 * circumradius(L)
  breaks <- handle.r.b.args(r, NULL, W, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max
  #
  if(correction == "Ang") {
    fname <- c("K", "list(L, I, J)")
    ylab <- quote(K[L,I,J](r))
  } else {
    fname <- c("K", "list(net, I, J)")
    ylab <- quote(K[net,I,J](r))
  }
  #
  if(np < 2) {
    # no pairs to count: return zero function
    zeroes <- rep(0, length(r))
    df <- data.frame(r = r, est = zeroes)
    K <- fv(df, "r", ylab,
            "est", . ~ r, c(0, rmax),
            c("r", makefvlabel(NULL, "hat", fname)),
            c("distance argument r", "estimated %s"),
            fname = fname)
    return(K)
  }
  #
  nI <- sum(I)
  nJ <- sum(J)
  whichI <- which(I)
  whichJ <- which(J)
  clash <- I & J
  has.clash <- any(clash)
  # compute pairwise distances
  if(exists("crossdist.lpp")) {
    DIJ <- crossdist(X[I], X[J], check=FALSE)
    if(has.clash) {
      # exclude pairs of identical points from consideration
      Iclash <- which(clash[I])
      Jclash <- which(clash[J])
      DIJ[cbind(Iclash,Jclash)] <- Inf
    }
  } else {
    D <- pairdist(X)
    diag(D) <- Inf
    DIJ <- D[I, J]
  }
  #---  compile into K function ---
  if(correction == "none" && is.null(reweight)) {
    # no weights (Okabe-Yamada)
    K <- compileK(DIJ, r, denom=denom, check=FALSE, fname=fname)
    K <- rebadge.as.crossfun(K, "K", "net", "I", "J")
    unitname(K) <- unitname(X)
    attr(K, "correction") <- correction
    return(K)
  }
  if(correction == "none")
     edgewt <- 1
  else {
     # inverse m weights (Ang's correction)
     # compute m[i,j]
     m <- matrix(1, nI, nJ)
     XPI <- XP[I]
     if(!has.clash) {
       for(k in seq_len(nJ)) {
         j <- whichJ[k]
         m[,k] <- countends(L, XPI, DIJ[, k])
       }
     } else {
       # don't count identical pairs
       for(k in seq_len(nJ)) {
         j <- whichJ[k]
         inotj <- (whichI != j)
         m[inotj, k] <- countends(L, XPI[inotj], DIJ[inotj, k])
       }
     }
     edgewt <- 1/m
  }
  # compute K
  wt <- if(!is.null(reweight)) edgewt * reweight else edgewt
  K <- compileK(DIJ, r, weights=wt, denom=denom, check=FALSE, fname=fname)
  ## rebadge and tweak
  K <- rebadge.as.crossfun(K, "K", "L", "I", "J")
  fname <- attr(K, "fname")
  # tack on theoretical value
  K <- bind.fv(K, data.frame(theo=r),
               makefvlabel(NULL, NULL, fname, "pois"),
               "theoretical Poisson %s")
  ## 
  unitname(K) <- unitname(X)
  fvnames(K, ".") <- rev(fvnames(K, "."))
  # show working
  if(showworking)
    attr(K, "working") <- list(DIJ=DIJ, wt=wt)
  attr(K, "correction") <- correction
  return(K)
}

