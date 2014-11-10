#
# exactMPLEstrauss.R
#
# 'exact' MPLE for stationary Strauss process
#
#  $Revision: 1.6 $  $Date: 2014/11/10 07:39:41 $
#

exactMPLEstrauss <- local({

  # main function
  exactMPLEstrauss <- function(X, R, ngrid=2048, plotit=FALSE, project=TRUE) {
#    n <- npoints(X)
    W <- as.owin(X)
    # border correction
    WminR <- erosion(W, R)
    bR <- (bdist.points(X) >= R)
    nR <- sum(bR)
    # evaluate neighbour counts for data points
    Tcounts <- crosspaircounts(X, X, R) - 1
    sumT  <- sum(Tcounts[bR])
    # determine the coefficients a_k for k = 0, 1, ...
    Z <- scanmeasure(X, R, dimyx=ngrid)
    Z <- Z[WminR, drop=FALSE]
    kcounts <- tabulate(as.vector(Z$v) + 1L)
    pixarea <- with(Z, xstep * ystep)
    A <- kcounts * pixarea
    # find optimal log(gamma)
    op <- optim(log(0.5), lpl, sco, method="L-BFGS-B",
                control=list(fnscale=-1),
                lower=-Inf, upper=if(project) 0 else Inf,
                A=A, sumT=sumT, nR=nR)
    loggamma <- op$par
    # plot?
    if(plotit) {
      x <- seq(log(1e-4), if(project) 0 else log(1e4), length=512)
      plot(x, lpl(x, A, sumT, nR),
           type="l",
           xlab=expression(log(gamma)),
           ylab=expression(log(PL(gamma))))
      abline(v=loggamma, lty=3)
    }
    # derive optimal beta 
    kmax <-length(A) - 1
    polypart <- A %*% exp(outer(0:kmax, loggamma))
    beta <- nR/polypart
    logbeta <- log(beta)
    result <- c(logbeta, loggamma)
    names(result) <- c("(Intercept)", "Interaction")
    return(result)
  }

  # helper functions (vectorised)
  # log pseudolikelihood
  lpl <- function(theta, A=A, sumT=sumT, nR=nR) {
    kmax <-length(A) - 1
    polypart <- A %*% exp(outer(0:kmax, theta))
    nR * (log(nR) - log(polypart) - 1) + theta * sumT
  }
  # pseudoscore
  sco <- function(theta, A=A, sumT=sumT, nR=nR) {
    kmax <- length(A) - 1
    kseq <- 0:kmax
    mat <- exp(outer(kseq, theta))
    polypart <- A %*% mat
    Dpolypart <- (A * kseq) %*% mat
    sumT - nR * Dpolypart/polypart
  }

  exactMPLEstrauss
})
