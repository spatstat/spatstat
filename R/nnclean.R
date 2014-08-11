#
#  nnclean.R
#
# Nearest-neighbour clutter removal
#
# Adapted from statlib file NNclean.q
# Authors: Simon Byers and Adrian Raftery
#
#  $Revision: 1.5 $   $Date: 2011/05/18 08:08:21 $
#

nnclean <- function(X, k, ...) {
  UseMethod("nnclean")
}

nnclean.pp3 <- function(X, k, ...,
                        convergence = 0.001, plothist = FALSE,
                        verbose=TRUE, maxit=50)
{
  # Adapted from statlib file NNclean.q
  # Authors: Simon Byers and Adrian Raftery
  # Adapted for spatstat by Adrian Baddeley

  stopifnot(inherits(X, "pp3"))
  validposint(k, "nnclean.pp3")

  kthNND <- nndist(X, k=k)  
  
  # apply classification algorithm
  em <- nncleanEngine(kthNND, k=k, d=3, ...,
                     tol=convergence, plothist=plothist,
                     verbose=verbose, maxit=maxit)

  # tack results onto point pattern as marks
  pp <- em$probs
  zz <- factor(em$z, levels=c(0,1))
  levels(zz) <- c("noise", "feature")
  mm <- hyperframe(prob=pp, label=zz)
  marks(X) <- cbind(marks(X), mm)
  return(X)
}

nnclean.ppp <-
  function(X, k, ...,
           edge.correct = FALSE, wrap = 0.1,
           convergence = 0.001, plothist = FALSE,
           verbose=TRUE, maxit=50)
{
  # Adapted from statlib file NNclean.q
  # Authors: Simon Byers and Adrian Raftery
  # Adapted for spatstat by Adrian Baddeley

  n <- X$n
  validposint(k, "nnclean.ppp")

  if(!edge.correct) {
    # compute vector of k-th nearest neighbour distances
    kthNND <- nndist(X, k=k)
  } else {
    # replicate data periodically
    # (ensuring original points are listed first)
    Xbox <- X[as.rectangle(X)]
    Xpand <- periodify(Xbox, ix=c(0,-1,1), iy=c(0,-1,1), check=FALSE)
    # trim to margin
    W <- expand.owin(X$window, (1+2*wrap)^2)
    Xpand <- Xpand[W]
    kthNND <- nndist(Xpand, k=k)
  }

  # apply classification algorithm
  em <- nncleanEngine(kthNND, k=k, d=2, ...,
                     tol=convergence, plothist=plothist,
                     verbose=verbose, maxit=maxit)

  # extract results
  pp <- em$probs
  zz <- em$z
  zz <- factor(zz, levels=c(0,1))
  levels(zz) <- c("noise", "feature")
  df <- data.frame(class=zz,prob=pp) 

  if(edge.correct) {
    # trim back to original point pattern
    df <- df[seq_len(X$n), ]
  }
  
  # tack on
  marx <- marks(X, dfok=TRUE)
  if(is.null(marx))
    marks(X, dfok=TRUE) <- df
  else 
    marks(X, dfok=TRUE) <- cbind(df, marx)

  return(X)
}

nncleanEngine <-
  function(kthNND, k, d, ..., 
           tol = 0.001, plothist = FALSE,
           verbose=TRUE, maxit=50)
{
  # Adapted from statlib file NNclean.q
  # Authors: Simon Byers and Adrian Raftery
  # Adapted for spatstat by Adrian Baddeley
  
  n <- length(kthNND)
  
  alpha.d <- (2. * pi^(d/2.))/(d * gamma(d/2.))

  # raise to power d for efficiency
  kNNDpowd <- kthNND^d
  
  #
  # Now use kthNND in E-M algorithm.
  # First set up starting guesses.
  #
  #
  probs <- rep(0, n)
  thresh <- (min(kthNND) + diff(range(kthNND))/3.)
  high <- (kthNND > thresh)
  delta <- ifelse(high, 1, 0)
  p <- 0.5
  lambda1 <- k/(alpha.d * mean(kNNDpowd[!high]))
  lambda2 <- k/(alpha.d * mean(kNNDpowd[ high]))
  loglik.old <- 0.
  loglik.new <- 1.
  #
  # Iterator starts here, 
  #
  Z <- !kthNND
  niter <- 0
  while(abs(loglik.new - loglik.old)/(1 + abs(loglik.new)) > tol) {
    if(niter >= maxit) {
      warning(paste("E-M algorithm failed to converge in",
                    maxit, "iterations"))
      break
    }
    niter <- niter + 1
    # E - step
    f1 <- dknn(kthNND[!Z], lambda=lambda1, k = k, d = d)
    f2 <- dknn(kthNND[!Z], lambda=lambda2, k = k, d = d)
    delta[!Z] <- deltaNZ <- (p * f1)/(p * f1 + (1 - p) * f2)
    delta[Z] <- 0
    # M - step
    p <- sum(delta)/n
    lambda1 <- (k * sum(delta))/(alpha.d * sum(kNNDpowd * delta))
    lambda2 <- (k * sum((1. - delta)))/(alpha.d * sum(kNNDpowd * (1. - delta)))
    # evaluate loglikelihood
    loglik.old <- loglik.new
    loglik.new <- sum( - p * lambda1 * alpha.d * (kNNDpowd * delta)
                      - (1. - p) * lambda2 * alpha.d * (kNNDpowd * (1 - delta))
                      + delta * k * log(lambda1 * alpha.d) +
			(1. - delta) * k * log(lambda2 * alpha.d))
    if(verbose) 
      cat(paste("Iteration", niter, "\tlogLik =", loglik.new,
                "\tp =", signif(p,4), "\n"))
  }
  if(plothist) {
    xlim <- c(0, max(kthNND))
    barheights <- hist(kthNND, nclass=40, plot=FALSE)$density
    support <- seq(from=xlim[1], to=xlim[2], length.out = 200.)
    fittedy <- p * dknn(support, lambda=lambda1, k = k, d = d) +
      (1. - p) * dknn(support, lambda=lambda2, k = k, d = d)
    ylim <- range(c(0, barheights, fittedy))
    xlab <- paste("Distance to", ordinal(k), "nearest neighbour")
    hist(kthNND,
         nclass=40,
         probability = TRUE,
         xlim = xlim, ylim=ylim, axes = TRUE,
         xlab = xlab, ylab = "Probability density")
    box()
    lines(support, fittedy, col="green")
  }
  #
  delta1 <- dknn(kthNND[!Z], lambda=lambda1, k = k, d = d)
  delta2 <- dknn(kthNND[!Z], lambda=lambda2, k = k, d = d)
  probs[!Z] <- delta1/(delta1 + delta2)
  probs[Z] <- 1
  #
  if(verbose) {
    cat("Estimated parameters:\n")
    cat(paste("p [cluster] =", signif(p, 5), "\n"))
    cat(paste("lambda [cluster] =", signif(lambda1, 5), "\n"))
    cat(paste("lambda [noise]   =", signif(lambda2, 5), "\n"))
  }
  #
  # z will be the classifications. 1= in cluster. 0= in noise. 
  #
  return(list(z = round(probs), probs = probs, lambda1 = lambda1, lambda2 = 
       lambda2, p = p, kthNND = kthNND, d=d, n=n, k=k))
}

