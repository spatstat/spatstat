#'
#'    sdr.R
#'
#'  Sufficient Dimension Reduction
#'
#'  Matlab original: Yongtao Guan
#'  Translated to R by:  Suman Rakshit
#'  Adapted for spatstat: Adrian Baddeley
#'
#'  GNU Public Licence 2.0 || 3.0
#'
#'    $Revision: 1.14 $  $Date: 2019/01/08 07:46:46 $
#'

sdr <- function(X, covariates, ...) {
  UseMethod("sdr")
}

sdr.ppp <- local({

  sdr.ppp <- function(X, covariates,
                      method=c("DR", "NNIR", "SAVE", "SIR", "TSE"),
                      Dim1=1, Dim2=1, predict=FALSE, ...) {
    stopifnot(is.ppp(X))
    method <- match.arg(method)
    trap.extra.arguments(...)
    #' ensure 'covariates' is a list of compatible images
    if(!inherits(covariates, "imlist") && !all(sapply(covariates, is.im)))
      stop("Argument 'covariates' must be a list of images")
    nc <- length(covariates)
    if(nc == 0)
      stop("Need at least one covariate!")
    if(nc < Dim1 + (method == "TSE") * Dim2)
      stop(paste(if(method == "TSE") "Dim1 + Dim2" else "Dim1",
                 "must not exceed the number of covariates"),
           call.=FALSE)
    if(nc > 1 && !do.call(compatible, unname(covariates)))
      covariates <- do.call(harmonise, covariates)
    #' extract corresponding pixel values including NA's
    Ypixval <- sapply(lapply(covariates, as.matrix), as.vector)
    #' compute sample mean and covariance matrix
    m <- colMeans(Ypixval, na.rm=TRUE)
    V <- cov(Ypixval, use="complete")
    #' evaluate each image at point data locations
    YX <- sapply(covariates, safelook, Y=X)
    #' apply precomputed standardisation
    Zx <- t(t(YX) - m) %*% matrixinvsqrt(V)
    #' ready
    coordsX <- coords(X)
    result <-
      switch(method,
             DR   =   calc.DR(COV=V, z=Zx,              Dim=Dim1),
             NNIR = calc.NNIR(COV=V, z=Zx, pos=coordsX, Dim=Dim1),
             SAVE = calc.SAVE(COV=V, z=Zx,              Dim=Dim1),
             SIR  =  calc.SIR(COV=V, z=Zx                      ),
             TSE  =  calc.TSE(COV=V, z=Zx, pos=coordsX, Dim1=Dim1, Dim2=Dim2)
             )
    #'
    covnames <- names(covariates) %orifnull% paste0("Y", 1:nc)
    dimnames(result$B) <- list(covnames, paste0("B", 1:ncol(result$B)))
    if(method == "TSE") {
      result$M1 <- namez(result$M1)
      result$M2 <- namez(result$M2)
    } else {
      result$M <- namez(result$M)
    }
    if(predict) result$Y <- sdrPredict(covariates, result$B)
    return(result)
  }

  safelook <- function(Z, Y, ...) { safelookup(Z, Y, ...) }

  namez <- function(M, prefix="Z") {
    dimnames(M) <- list(paste0(prefix, 1:nrow(M)),
                     paste0(prefix, 1:ncol(M)))
    return(M)
  }
  
  sdr.ppp
})

sdrPredict <- function(covariates, B) {
  if(!is.matrix(B)) {
    if(is.list(B) && is.matrix(BB <- B$B)) B <- BB else
    stop("B should be a matrix, or the result of a call to sdr()",
         call.=FALSE)
  }
  if(!inherits(covariates, "imlist") && !all(sapply(covariates, is.im)))
    stop("Argument 'covariates' must be a list of images")
  stopifnot(nrow(B) == length(covariates))
  result <- vector(mode="list", length=ncol(B))
  for(j in seq_along(result)) {
    cj <- as.list(B[,j])
    Zj <- mapply("*", cj, covariates, SIMPLIFY=FALSE)
    result[[j]] <- im.apply(Zj, sum)
  }
  names(result) <- colnames(B)
  return(as.solist(result))
}

##............ DR (Directional Regression) ..........................
calc.DR <- function(COV, z, Dim){
  ## Description: Naive Directional Regression Method
  ## Input:
  ##   COV - cov{X(s)}
  ##   z   - standardized X(s) on SPP locations
  ##   Dim - the CS dimension
  ## Output:
  ##   B   - the estimated CS basis
  ##   M - the kernel matrix
  ss <- nrow(z)
  ncov <- ncol(z)
  M1 <- (t(z) %*% z)/ss - diag(1,ncov)
  M1 <- M1 %*% M1                             # the SAVE kernel
  covMean <- matrix(colMeans(z),ncol=1)
  M2 <- covMean %*% t(covMean)
  M3 <- M2 * (base::norm(covMean, type="2"))^2      # the SIR kernel
  M2 <- M2 %*% M2                             # the SIR-2 kernel
  M <- (M1 + M2 + M3)/3                       # the DR kernel
  SVD <- svd(M)
  B <- SVD$u[,1:Dim]
  B <- matrixinvsqrt(COV) %*% B               # back to original scale
  return(list(B=B, M=M))
}


## ............ NNIR (Nearest Neighbor Inverse Regression) ...........
calc.NNIR <- function(COV, z, pos, Dim) {
  ## Description: Nearest Neighbor Inverse Regression
  ## Input:
  ##   COV - cov{X(s)}
  ##   z   - standardized X(s) on SPP locations
  ##   pos - the position of SPP events
  ##   Dim - the CS dimension
  ## Output:
  ##   B   - the estimated CS basis
  ##   M - the kernel matrix

  ss   <- nrow(z)   # sample size
#  ncov <- ncol(z)   # predictor dimension
  jj <- nnwhich(pos) # identify nearest neighbour of each point
  dir <- z - z[jj, , drop=FALSE] # empirical direction
  IM <- sumouter(dir) # inverse of kernel matrix: sum of outer(dir[i,], dir[i,])
  M <- solve(IM/ss)   # invert kernel matrix
  SVD <- svd(M)
  B <- matrixinvsqrt(COV) %*% SVD$u[, 1:Dim, drop=FALSE]
  return(list(B=B, M=M))
}

## ...........  SAVE (Sliced Average Variance Estimation) ...........
calc.SAVE <- function(COV, z, Dim){
  ## Description: Naive Directional Regression Method
  ## Input
  ##   COV - cov{X(s)}
  ##   z - standardized X(s) on SPP locations
  ##   Dim - the central space dimension
  ## Value
  ##   B - the estimated CS basis
  ##   M - the kernel matrix
#  ss <- nrow(z)
  ncov <- ncol(z)
  M <- diag(1,ncov) - cov(z)
  M <- M %*% M
  SVD <- svd(M)
  B <- SVD$u[,1:Dim]
  B <- matrixinvsqrt(COV) %*% B
  return(list(B=B, M=M))
}

##..........  SIR (Sliced Inverse Regression) ......................
calc.SIR <- function(COV, z){
  ## Description: Naive Directional Regression Method
  ## Input:
  ##   COV - cov{X(s)}
  ##   z   - standardized X(s) on SPP locations
  ## Output:
  ##   B   - the estimated CS basis
  ##   M - the kernel matrix
  covMean <- colMeans(z)
  B <- matrixinvsqrt(COV) %*% covMean    # do SIR estimation
  B <- B/sqrt(sum(B^2))                   # normalise to unit length
  M <- covMean %*% t(covMean)             # create kernel matrix
  return(list(B=B, M=M))
}

## .............  TSE (Two-Step Estimation) ....................
calc.TSE <- function(COV, z, pos, Dim1, Dim2) {
  ## Description: A Two-Step Method
  ## Input:
  ##   COV - cov{X(s)}
  ##   z   - standardized X(s) on SPP locations
  ##   Dim1 - the S1 dimension
  ##   Dim2 - the S2 dimension
  ## Output:
  ##   B   - the estimated CS basis. Its first Dim1 columns
  ##       are estimating S1 and the remaining Dim2 columns are
  ##       estimating S2. In case of null space, a zero vector is reported.
  ##   M1  - the kernel matrix of DR
  ##   M2  - the kernel matrix of NNIR, which might be subject 
  ##         to some change, depending on the results of M1.

#  ss   <- nrow(z)  # sample size
  ncov <- ncol(z)  # predictor dimension

  est1 <- calc.DR(COV, z, ncov)           # do DR estimation
  est2 <- calc.NNIR(COV, z, pos, ncov)  # do NNIR estimation
  M1 <- est1$M
  M2 <- est2$M

  if(Dim1 > 0) {
    U <- svd(M1)$u
    B1 <- U[ , 1:Dim1, drop=FALSE]  # get S1 estimate
    Q  <- diag(1, ncov) - B1 %*% solve(t(B1) %*% B1) %*% t(B1)
                     # contract orthogonal basis
    M2 <- Q %*% M2 %*% Q  # do constrained NNIR
  } else {
    B1 <- matrix(0, ncov, 1)
  }

  if(Dim2 > 0) {
    U <- svd(M2)$u   # do SVD for possibly updated M2
    B2 <- U[ , 1:Dim2, drop=FALSE]  # get basis estimator
  } else {
    B2 <- matrix(0, ncov, 1)
  }
  B <- matrixinvsqrt(COV) %*% cbind(B1,B2)
  return(list(B=B, M1=M1, M2=M2))
}

## //////////////////  ADDITIONAL FUNCTIONS /////////////////////

subspaceDistance <- function(B0,B1) {
  ## ======================================================== #
  ## Evaluate the distance between the two linear spaces S(B0) and S(B1). 
  ## The measure used is the one proposed by Li et al. (2004). 
  ## ======================================================== #
  stopifnot(is.matrix(B0))
  stopifnot(is.matrix(B1))
  Proj0 <- B0 %*% solve((t(B0) %*% B0)) %*% t(B0)  # Proj matrix on S(B0)
  lam <- svd(B1) # check whether B1 is singular
  U <- lam$u
  D <- lam$d
#  V <- lam$v
  B2 <- U[, D > 1e-09] # keep non-singular directions
  Proj1 <- B2 %*% solve((t(B2) %*% B2)) %*% t(B2) # Proj matrix on S(B.hat)
  Svd <- svd(Proj0 - Proj1)  # Do svd for P0-P1
  dist <- max(abs(Svd$d)) # Get the maximum absolute svd value
  return(dist)
}

dimhat <- function(M){
  #' Description: Maximum Descent Estimator for CS Dim
  #' Input:
  #'   M   - the estimated kernel matrix
  #' Output:
  #'   dimhat   - the estimated CS dim (assume dim>0)
  stopifnot(is.matrix(M))
  ncov <- ncol(M)                       # predictor dimension
  maxdim <- max((ncov-1), 5)            # maximum structure dimension
  SVD <- svd(M)                         # svd of kernel matrix
  lam <- SVD$d
  eps <- 1e-06
  lam <- lam + rep(eps,ncov)            # add ridge effect
  lam1 <- lam[-ncov]
  lam2 <- lam[-1]
  dif <- lam1/lam2
  dif <- dif[1 : maxdim]                # the magnitude of drop
  retval <- which.max(dif)              # find Maximum Descent estimator
  return(retval)
}
