#
# linearKmulti
#
# $Revision: 1.5 $ $Date: 2013/01/31 02:46:00 $
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
  iname <- make.parseable(paste(i))
  result <-
    rebadge.fv(result,
               substitute(linearK[i ~ dot](r), list(i=iname)),
               paste("linearK[", iname, "~ symbol(\"\\267\")]"),
               new.yexp=substitute(linearK[i ~ symbol("\267")](r), list(i=iname)))
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
  iname <- make.parseable(paste(i))
  jname <- make.parseable(paste(j))
  result <-
    rebadge.fv(result, 
               substitute(linearKcross[i,j](r), list(i=iname,j=jname)),
               sprintf("linearK[list(%s,%s)]", iname, jname),
               new.yexp=substitute(linearK[list(i,j)](r),
                                   list(i=iname,j=jname)))
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
  switch(correction,
         Ang  = {
           ylab <- quote(Kmulti[L](r))
           fname <- "Kmulti[L]"
         },
         none = {
           ylab <- quote(Kmulti[net](r))
           fname <- "Kmulti[net]"
         })
  K <- rebadge.fv(K, new.ylab=ylab, new.fname=fname)
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
  iname <- make.parseable(paste(i))
  result <-
    rebadge.fv(result,
               substitute(linearK[inhom, i ~ dot](r), list(i=iname)),
               paste("linearK[list(inhom,", iname, "~ symbol(\"\\267\"))]"),
               new.yexp=substitute(linearK[list(inhom, i ~ symbol("\267"))](r), list(i=iname)))
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
  iname <- make.parseable(paste(i))
  jname <- make.parseable(paste(j))
  result <-
    rebadge.fv(result, 
               substitute(linearK[inhom,i,j](r), list(i=iname,j=jname)),
               sprintf("linearK[list(inhom,%s,%s)]", iname, jname),
               new.yexp=substitute(linearK[list(inhom,i,j)](r),
                                   list(i=iname,j=jname)))
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
  switch(correction,
         Ang  = {
           ylab <- quote(Kmulti[L](r))
           fname <- "Kmulti[L]"
         },
         none = {
           ylab <- quote(Kmulti[net](r))
           fname <- "Kmulti[net]"
         })
  K <- rebadge.fv(K, new.ylab=ylab, new.fname=fname)
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
  if(np < 2) {
    # no pairs to count: return zero function
    zeroes <- rep(0, length(r))
    df <- data.frame(r = r, est = zeroes)
    K <- fv(df, "r", substitute(linearK(r), NULL),
            "est", . ~ r, c(0, rmax),
            c("r", "%s(r)"),
            c("distance argument r", "estimated %s"),
            fname = "linearK")
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
    K <- compileK(DIJ, r, denom=denom, check=FALSE)
    unitname(K) <- unitname(X)
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
  K <- compileK(DIJ, r, weights=wt, denom=denom, check=FALSE)
  # tack on theoretical value
  K <- bind.fv(K, data.frame(theo=r), "%s[theo](r)",
               "theoretical Poisson %s")
  unitname(K) <- unitname(X)
  fvnames(K, ".") <- rev(fvnames(K, "."))
  # show working
  if(showworking)
    attr(K, "working") <- list(DIJ=DIJ, wt=wt)
  return(K)
}

