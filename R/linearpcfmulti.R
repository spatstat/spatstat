#
# linearpcfmulti.R
#
# $Revision: 1.3 $ $Date: 2014/10/24 00:22:30 $
#
# pair correlation functions for multitype point pattern on linear network
#
#

linearpcfdot <- function(X, i, r=NULL, ..., correction="Ang") {
  if(!is.multitype(X, dfok=FALSE)) 
	stop("Point pattern must be multitype")
  marx <- marks(X)
  lev <- levels(marx)
  if(missing(i) || is.null(i)) i <- lev[1] else
    if(!(i %in% lev)) stop(paste("i = ", i , "is not a valid mark"))  
  I <- (marx == i)
  J <- rep(TRUE, npoints(X))  # i.e. all points
  result <- linearpcfmulti(X, I, J,
                           r=r, correction=correction, ...)
  correction <- attr(result, "correction")
  type <- if(correction == "Ang") "L" else "net"
  result <- rebadge.as.dotfun(result, "g", type, i)
  return(result)
}

linearpcfcross <- function(X, i, j, r=NULL, ..., correction="Ang") {
  if(!is.multitype(X, dfok=FALSE)) 
	stop("Point pattern must be multitype")
  marx <- marks(X)
  lev <- levels(marx)
  if(missing(i) || is.null(i)) i <- lev[1] else
    if(!(i %in% lev)) stop(paste("i = ", i , "is not a valid mark"))
  if(missing(j) || is.null(j)) j <- lev[2] else
    if(!(j %in% lev)) stop(paste("j = ", j , "is not a valid mark"))
  #
  if(i == j) {
    result <- linearpcf(X[marx == i], r=r, correction=correction, ...)
  } else {
    I <- (marx == i)
    J <- (marx == j)
    result <- linearpcfmulti(X, I, J, r=r, correction=correction, ...)
  }
  # rebrand
  correction <- attr(result, "correction")
  type <- if(correction == "Ang") "L" else "net"
  result <- rebadge.as.crossfun(result, "g", type, i, j)
  return(result)
}

linearpcfmulti <- function(X, I, J, r=NULL, ..., correction="Ang") {
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
  # compute pcf
  denom <- (nI * nJ - nIandJ)/lengthL
  g <- linearPCFmultiEngine(X, I, J, r=r, denom=denom, correction=correction, ...)
  # set appropriate y axis label
  correction <- attr(g, "correction")
  type <- if(correction == "Ang") "L" else "net"
  g <- rebadge.as.crossfun(g, "g", type, "I", "J")
  attr(g, "correction") <- correction
  return(g)
}

# ................ inhomogeneous ............................

linearpcfdot.inhom <- function(X, i, lambdaI, lambdadot,
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
  result <- linearpcfmulti.inhom(X, I, J, lambdaI, lambdadot, 
                               r=r, correction=correction, normalise=normalise,
                               ...)
  correction <- attr(result, "correction")
  type <- if(correction == "Ang") "L, inhom" else "net, inhom"
  result <- rebadge.as.dotfun(result, "g", type, i)
  return(result)
}

linearpcfcross.inhom <- function(X, i, j, lambdaI, lambdaJ,
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
    result <- linearpcfinhom(X[I], lambda=lambdaI, r=r,
                           correction=correction, normalise=normalise, ...)
  } else {
    I <- (marx == i)
    J <- (marx == j)
    result <- linearpcfmulti.inhom(X, I, J, lambdaI, lambdaJ,
                                 r=r, correction=correction,
                                 normalise=normalise, ...)
  }
  # rebrand
  correction <- attr(result, "correction")
  type <- if(correction == "Ang") "L, inhom" else "net, inhom"
  result <- rebadge.as.crossfun(result, "g", type, i, j)
  return(result)
}

linearpcfmulti.inhom <- function(X, I, J, lambdaI, lambdaJ,
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

  # compute pcf
  weightsIJ <- outer(1/lambdaI, 1/lambdaJ, "*")
  denom <- if(!normalise) lengthL else sum(1/lambdaI)
  g <- linearPCFmultiEngine(X, I, J, r=r,
                            reweight=weightsIJ, denom=denom,
                            correction=correction, ...)
  # set appropriate y axis label
  correction <- attr(g, "correction")
  type <- if(correction == "Ang") "L, inhom" else "net, inhom"
  g <- rebadge.as.crossfun(g, "g", type, "I", "J")
  attr(g, "correction") <- correction
  return(g)
}

# .............. internal ...............................

linearPCFmultiEngine <- function(X, I, J, ..., r=NULL, reweight=NULL, denom=1,
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
    fname <- c("g", "list(L, I, J)")
    ylab <- quote(g[L,I,J](r))
  } else {
    fname <- c("g", "list(net, I, J)")
    ylab <- quote(g[net,I,J](r))
  }
  #
   if(np < 2) {
    # no pairs to count: return zero function
    zeroes <- rep(0, length(r))
    df <- data.frame(r = r, est = zeroes)
    g <- fv(df, "r", ylab,
            "est", . ~ r, c(0, rmax),
            c("r", makefvlabel(NULL, "hat", fname)), 
            c("distance argument r", "estimated %s"),
            fname = fname)
    unitname(g) <- unitname(X)
    attr(g, "correction") <- correction
    return(g)
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
  #---  compile into pair correlation function ---
  if(correction == "none" && is.null(reweight)) {
    # no weights (Okabe-Yamada)
    g <- compilepcf(DIJ, r, denom=denom, check=FALSE, fname=fname)
    g <- rebadge.as.crossfun(g, "g", "net", "I", "J")    
    unitname(g) <- unitname(X)
    attr(g, "correction") <- correction
    return(g)
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
  # compute pcf
  wt <- if(!is.null(reweight)) edgewt * reweight else edgewt
  g <- compilepcf(DIJ, r, weights=wt, denom=denom, check=FALSE, ...,
                  fname=fname)
  ## rebadge and tweak
  g <- rebadge.as.crossfun(g, "g", "L", "I", "J")
  fname <- attr(g, "fname")
  # tack on theoretical value
  g <- bind.fv(g, data.frame(theo=rep(1,length(r))),
               makefvlabel(NULL, NULL, fname, "pois"),
               "theoretical Poisson %s")
  unitname(g) <- unitname(X)
  fvnames(g, ".") <- rev(fvnames(g, "."))
  # show working
  if(showworking)
    attr(g, "working") <- list(DIJ=DIJ, wt=wt)
  attr(g, "correction") <- correction
  return(g)
}

