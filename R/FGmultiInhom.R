#'
#'     FGmultiInhom.R
#' 
#'     inhomogeneous multitype G and F functions
#'
#'     Original code by Ottmar Cronie and Marie-Colette van Lieshout
#'
#'     Rewritten for spatstat by Adrian Baddeley
#'
#'     GmultiInhom
#'     FmultiInhom
#'
#'      $Revision: 1.5 $ $Date: 2016/07/06 07:01:22 $

GmultiInhom <- function(X, I, J, 
                        lambda=NULL, lambdaI=NULL, lambdaJ=NULL,
                        lambdamin=NULL,
                        ...,
                        r=NULL, 
                        ReferenceMeasureMarkSetI=NULL,
                        ratio=FALSE){
  if(!is.ppp(X) || !is.marked(X))
    stop("X should be a marked point pattern")
  W <- Window(X)
  nX <- npoints(X)
  
  #' handle r argument
  rmax <- rmax.rule("G", W, intensity(X))
  bks <- handle.r.b.args(r, NULL, W, rmaxdefault=rmax)
  r    <- bks$r
  rmax <- bks$max
  nr   <- length(r)
  
  #' Accept any kind of index for I; convert it to a logical index
  I <- ppsubset(X, I)
  if(is.null(I))
    stop("I must be a valid subset index")
  XI <- X[I]
  nI <- sum(I)
  if (nI == 0) 
    stop("No points satisfy condition I")
  
  if(!is.null(ReferenceMeasureMarkSetI)) {
    check.1.real(ReferenceMeasureMarkSetI)
    stopifnot(ReferenceMeasureMarkSetI >= 0)
  }

  #' likewise for J
  if(missing(J) || is.null(J)) {
    J <- rep(TRUE, nX)
  } else {
    J <- ppsubset(X, J)
  }
  XJ <- X[J]
  nJ <- sum(J)
  if (nJ == 0) 
    stop("No points satisfy condition J")

  #' supply either lambda, or lambdaI and lambdaJ
  lam.given <- !is.null(lambda)
  lamIJ.given <- !is.null(lambdaI) || !is.null(lambdaJ)
  if(lam.given == lamIJ.given || is.null(lambdaI) != is.null(lambdaJ))
    stop(paste("Supply either a vector lambda of length equal to npoints(X),",
               "or two vectors lambdaI, lambdaJ of lengths",
               "equal to npoints(X[I]) and npoints(X[J]) respectively"),
         call.=FALSE)
  
  if(lamIJ.given) {
    #' lambdaI and lambdaJ given
    check.nvector(lambdaI, nI, things="points of X[I]")
    stopifnot(all(lambdaI > 0))
    check.nvector(lambdaJ, nJ, things="points of X[J]")
    stopifnot(all(lambdaJ > 0))
    if(is.null(lambdamin)){
      stop(paste("Supply lambdamin - a single positive number which is",
                 "smaller than the values in lambdaJ"),
           call.=FALSE)
    }
    check.1.real(lambdamin)
    stopifnot(lambdamin > 0)
    stopifnot(lambdamin <= min(lambdaJ))
  } else {
    #' lambda given
    check.nvector(lambda, nX, things="points of X")
    stopifnot(all(lambda > 0))
    lambdaI <- lambda[I]
    lambdaJ <- lambda[J]
    if(is.null(lambdamin)){
      stop(paste("Supply lambdamin - a single positive number which is",
                 "smaller than the values in lambda"),
           call.=FALSE)
    }
    check.1.real(lambdamin)
    stopifnot(lambdamin > 0)
    stopifnot(lambdamin <= min(lambda))
  }
  
  #' Calculate 1/lambda(x_i,y_i,m_i))
  #'           for all (x_i,y_i,m_i) with m_i in I
  invlambdaI <- 1/lambdaI
  #' Calculate (1 - lambda_min/lambda(x_i,y_i,m_i))
  #'           for all (x_i,y_i,m_i) with m_i in J
  Coeff <- 1-(lambdamin/lambdaJ)
  ## CoeffMatrix <- matrix(rep(Coeff,times=nI), nrow=nI, byrow=TRUE)

  #' distances
  ## DistanceXItoXJ <- crossdist(XI,XJ)

  #' eroded areas and boundary distances
  areaWr <- eroded.areas(W, r)
  bdistXI <- bdist.points(XI)

  #' for each point x in XI, determine largest r such that x \in W-r
  ibI <- fastFindInterval(bdistXI, r, labels=TRUE)
  #' count of points inside W-r for each r
  ## NumberEroded <- revcumsum(table(ibI))
    
  #' denominator
  #' sum invlambdaI for all points x \in W-r
  DenominatorN <- c(sum(invlambdaI),
                    revcumsum(natozero(tapply(invlambdaI, ibI, sum))))
  if(!is.null(ReferenceMeasureMarkSetI))
    DenominatorA <- areaWr * ReferenceMeasureMarkSetI

  #' local products of weights
  #' sort data points in order of increasing x coordinate
  xxI <- XI$x
  yyI <- XI$y
  oXI <- fave.order(xxI)
  xIord <- xxI[oXI]
  yIord <- yyI[oXI]
  #'
  xxJ <- XJ$x
  yyJ <- XJ$y
  vvJ <- Coeff
  oXJ <- fave.order(xxJ)
  xJord <- xxJ[oXJ]
  yJord <- yyJ[oXJ]
  vJord <- vvJ[oXJ]
  # compute local cumulative products
  z <- .C("locxprod",
          ntest = as.integer(nI),
          xtest = as.double(xIord),
          ytest = as.double(yIord),
          ndata = as.integer(nJ),
          xdata = as.double(xJord),
          ydata = as.double(yJord),
          vdata = as.double(vJord),
          nr = as.integer(nr),
          rmax = as.double(rmax),
          ans = as.double(numeric(nI * nr)))
  ans <- matrix(z$ans, nrow=nr, ncol=nI)
  #' revert to original ordering
  loccumprod <- matrix(,  nrow=nr, ncol=nI)
  loccumprod[, oXI] <- ans

  #' border correction
  outside <- outer(r, bdistXI, ">")
  loccumprod[outside] <- 0
  #' weight by 1/lambdaI
  wlcp <- loccumprod * matrix(invlambdaI, byrow=TRUE, nr, nI)
  #' sum over I for each fixed r
  numer <- .rowSums(wlcp, nr, nI)

  # pack up
  Gdf <- data.frame(r=r, theo = 1 - exp(- lambdamin * pi * r^2))
  desc <- c("distance argument r", "theoretical Poisson %s")
  theo.denom <- rep.int(nI, nr)
  fname <- c("G", "list(inhom,I,J)")
  G <- ratfv(Gdf, NULL, theo.denom,
             "r", quote(G[inhom, I, J](r)),
             "theo", NULL, c(0,rmax),
             c("r", makefvlabel(NULL, NULL, fname, "pois")),
             desc,
             fname=fname,
             yexp=quote(G[list(inhom,I,J)](r)),
             ratio=ratio)
  # add border corrected (Hamilton Principle) estimate
  G <- bind.ratfv(G,
                  data.frame(bord=DenominatorN-numer), DenominatorN,
                  makefvlabel(NULL, "hat", fname, "bord"),
                  "border estimate of %s",
                  "bord",
                  ratio=ratio)
  fvnames(G, ".") <- c("bord", "theo")
  # add modified border corrected (non-Hamilton-Principle) estimate
  if(!is.null(ReferenceMeasureMarkSetI)) {
    G <- bind.ratfv(G,
                    data.frame(bordm=DenominatorA-numer),
                    DenominatorA,
                    makefvlabel(NULL, "hat", fname, "bordm"),
                    "modified border estimate of %s",
                    "bordm",
                    ratio=ratio)
    fvnames(G, ".") <- c("bord", "bordm", "theo")
  }
  # 
  formula(G) <- . ~ r
  unitname(G) <- unitname(X)
  if(ratio)
    G <- conform.ratfv(G)

  return(G)
}

#' marked inhomogeneous F

FmultiInhom <- function(X, J,
                        lambda=NULL,lambdaJ=NULL,
                        lambdamin=NULL,
                        ...,
                        r=NULL) {
  if(!is.ppp(X) || !is.marked(X))
    stop("X should be a marked point pattern")
  nX <- npoints(X)
  
  #' Accept any kind of index for J; convert it to a logical index
  J <- ppsubset(X, J)
  if(is.null(J))
    stop("J must be a valid subset index")
  XJ <- X[J]
  nJ <- sum(J)
  if (nJ == 0) 
    stop("No points satisfy condition J")
  
  if(is.null(lambda) == is.null(lambdaJ))
    stop(paste("Supply either a vector lambda of length equal to npoints(X),",
               "or a vector lambdaJ of length equal to npoints(X[J])"),
         call.=FALSE)
  if(is.null(lambdamin))
    stop("Supply a value for lambdamin", call.=FALSE)
  check.1.real(lambdamin)
  
  if(!is.null(lambda)) {
    check.nvector(lambda, nX)
    stopifnot(all(lambda > 0))
    stopifnot(lambdamin <= min(lambda[J]))
    lambdaJ <- lambda[J]
  } else {
    check.nvector(lambdaJ, nJ)
    stopifnot(all(lambdaJ > 0))
    stopifnot(lambdamin <= min(lambdaJ))
  }

  FJ <- Finhom(XJ, lambda=lambdaJ, lmin=lambdamin, r=r)
  FJ <- rebadge.fv(FJ,
                   new.ylab  = quote(F[inhom, J](r)),
                   new.fname = c("F", "list(inhom,J)"),
                   new.yexp   = quote(F[list(inhom,J)](r)))
  return(FJ)
}
