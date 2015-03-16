#
#  nnclean.R
#
# Nearest-neighbour clutter removal
#
# Adapted from statlib file NNclean.q
# Authors: Simon Byers and Adrian Raftery
#
#  $Revision: 1.14 $   $Date: 2014/11/10 10:55:29 $
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

  Xname <- short.deparse(substitute(X))
  
  stopifnot(inherits(X, "pp3"))
  validposint(k, "nnclean.pp3")

  kthNND <- nndist(X, k=k)  
  
  # apply classification algorithm
  em <- do.call(nncleanEngine,
                resolve.defaults(list(kthNND, k=k),
                                 list(...),
                                 list(d=3, tol=convergence, plothist=plothist,
                                      verbose=verbose, maxit=maxit,
                                      Xname=Xname)))

  # tack results onto point pattern as marks
  pp <- em$probs
  zz <- factor(em$z, levels=c(0,1))
  levels(zz) <- c("noise", "feature")
  mm <- hyperframe(prob=pp, label=zz)
  marks(X) <- cbind(marks(X), mm)
  attr(X, "theta") <- em[c("lambda1", "lambda2", "p")]
  attr(X, "info") <- em[c("d", "niter", "maxit", "converged")]
  attr(X, "hist") <- em$hist
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

  Xname <- short.deparse(substitute(X))
  
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
  em <- do.call(nncleanEngine,
                resolve.defaults(list(kthNND, k=k),
                                 list(...),
                                 list(d=2, tol=convergence, plothist=plothist,
                                      verbose=verbose, maxit=maxit,
                                      Xname=Xname)))

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

  attr(X, "theta") <- em[c("lambda1", "lambda2", "p")]
  attr(X, "info") <- em[c("d", "niter", "maxit", "converged")]
  attr(X, "hist") <- em$hist
  return(X)
}

nncleanEngine <-
  function(kthNND, k, d, ..., 
           tol = 0.001, maxit = 50,
           plothist = FALSE, lineargs = list(), 
           verbose=TRUE, Xname="X")
{
  ## Adapted from statlib file NNclean.q
  ## Authors: Simon Byers and Adrian Raftery
  ## Adapted for spatstat by Adrian Baddeley
  
  n <- length(kthNND)

  ## Undocumented extension by Adrian Baddeley 2014
  ## Allow different dimensions in feature and noise.
  ## d[1] is cluster dimension.
  
  d <- ensure2vector(d)
  alpha.d <- (2. * pi^(d/2.))/(d * gamma(d/2.))

  # raise to power d for efficiency
  kNNDpowd1 <- kthNND^(d[1])
  kNNDpowd2 <- kthNND^(d[2])
  
  #
  # Now use kthNND in E-M algorithm
  # First set up starting guesses.
  #
  #
  probs <- numeric(n)
  thresh <- (min(kthNND) + diff(range(kthNND))/3.)
  high <- (kthNND > thresh)
  delta <- as.integer(high)
  p <- 0.5
  lambda1 <- k/(alpha.d[1] * mean(kNNDpowd1[!high]))
  lambda2 <- k/(alpha.d[2] * mean(kNNDpowd2[ high]))
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
                    maxit, ngettext(maxit, "iteration", "iterations")),
              call.=FALSE)
      break
    }
    niter <- niter + 1
    # E - step
    f1 <- dknn(kthNND[!Z], lambda=lambda1, k = k, d = d[1])
    f2 <- dknn(kthNND[!Z], lambda=lambda2, k = k, d = d[2])
    delta[!Z] <- (p * f1)/(p * f1 + (1 - p) * f2)
    delta[Z] <- 0
    # M - step
    sumdelta <- sum(delta)
    negdelta <- 1. - delta
    p <- sumdelta/n
    lambda1 <- (k * sumdelta)/(alpha.d[1] * sum(kNNDpowd1 * delta))
    lambda2 <- (k * (n - sumdelta))/(alpha.d[2] * sum(kNNDpowd2 * negdelta))
    # evaluate marginal loglikelihood
    loglik.old <- loglik.new
    loglik.new <- sum( - p * lambda1 * alpha.d[1] * (kNNDpowd1 * delta)
                      - (1. - p) * lambda2 * alpha.d[2] * (kNNDpowd2 * negdelta)
                      + delta * k * log(lambda1 * alpha.d[1]) +
			negdelta * k * log(lambda2 * alpha.d[2]))
    if(verbose) 
      cat(paste("Iteration", niter, "\tlogLik =", loglik.new,
                "\tp =", signif(p,4), "\n"))
  }
  if(plothist) {
    dotargs <- list(...)
    if(spatstat.options('monochrome'))
      dotargs <- col.args.to.grey(dotargs)
    ## compute plot limits to include both histogram and density
    xlim <- c(0, max(kthNND))
    H <- do.call("hist",
                 resolve.defaults(list(kthNND, plot=FALSE, warn.unused=FALSE),
                                  dotargs,
                                  list(nclass=40)))
    barheights <- H$density
    support <- seq(from=xlim[1], to=xlim[2], length.out = 200)
    fittedy <- p * dknn(support, lambda=lambda1, k = k, d = d[1]) +
      (1 - p) * dknn(support, lambda=lambda2, k = k, d = d[2])
    ylim <- range(c(0, barheights, fittedy))
    xlab <- paste("Distance to", ordinal(k), "nearest neighbour")
    ## now plot it (unless overridden by plot=FALSE)
    reallyplot <- resolve.1.default("plot", list(...), list(plot=TRUE))
    H <- do.call("hist",
                 resolve.defaults(list(kthNND, probability=TRUE),
                                  dotargs,
                                  list(plot=TRUE,
                                       warn.unused=reallyplot,
                                       nclass=40,
                                       xlim = xlim, ylim=ylim,
                                       xlab = xlab,
                                       ylab = "Probability density",
                                       axes = TRUE, main="")))
    H$xname <- xlab
    if(reallyplot) {
      box()
      lineargs <- resolve.defaults(lineargs, list(col="green", lwd=2))
      if(spatstat.options("monochrome"))
        lineargs <- col.args.to.grey(lineargs)
      do.call("lines", append(list(x=support, y=fittedy), lineargs))
    }
  }
  #
  delta1 <- dknn(kthNND[!Z], lambda=lambda1, k = k, d = d[1])
  delta2 <- dknn(kthNND[!Z], lambda=lambda2, k = k, d = d[2])
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
  return(list(z = round(probs),
              probs = probs,
              lambda1 = lambda1, lambda2 = lambda2, p = p,
              kthNND = kthNND, d=d, n=n, k=k,
              niter = niter, maxit = maxit,
              converged = (niter >= maxit),
              hist=if(plothist) H else NULL))
}

