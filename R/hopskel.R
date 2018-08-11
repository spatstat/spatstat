##
## hopskel.R
##     Hopkins-Skellam test
##
##  $Revision: 1.2 $  $Date: 2014/09/23 08:24:36 $

hopskel <- function(X) {
  stopifnot(is.ppp(X))
  n <- npoints(X)
  if(n < 2) return(NA)
  dX <- nndist(X)
  U <- runifpoint(n, Window(X))
  dU <- nncross(U, X, what="dist")
  A <- mean(dX^2)/mean(dU^2)
  return(A)
}

hopskel.test <- function(X, ..., 
                         alternative=c("two.sided", "less", "greater",
                           "clustered", "regular"),
                         method=c("asymptotic", "MonteCarlo"),
                         nsim=999
                         ) {
  Xname <- short.deparse(substitute(X))

  verifyclass(X, "ppp")
  W <- Window(X)
  n <- npoints(X)

  method <- match.arg(method)
  
  # alternative hypothesis
  alternative <- match.arg(alternative)
  if(alternative == "clustered") alternative <- "less"
  if(alternative == "regular") alternative <- "greater"
  altblurb <-
    switch(alternative,
           two.sided="two-sided",
           less="clustered (A < 1)",
           greater="regular (A > 1)")

  ## compute observed value
  statistic <- hopskel(X)
  ## p-value
  switch(method,
         asymptotic = {
           ## F-distribution
           nn <- 2 * n
           p.value <-
             switch(alternative,
                    less = pf(statistic, nn, nn, lower.tail=TRUE),
                    greater = pf(statistic, nn, nn, lower.tail=FALSE),
                    two.sided = 2 *
                    pf(statistic, nn, nn, lower.tail=(statistic < 1)))
           pvblurb <- "using F distribution"
         },
         MonteCarlo = {
           ## Monte Carlo p-value
           sims <- numeric(nsim)
           for(i in 1:nsim) {
             Xsim <- runifpoint(n, win=W)
             sims[i] <- hopskel(Xsim)
             p.upper <- (1 + sum(sims >= statistic))/(1 + nsim)
             p.lower <- (1 + sum(sims <= statistic))/(1 + nsim)
             p.value <- switch(alternative,
                               less=p.lower,
                               greater=p.upper,
                               two.sided=2*min(p.lower, p.upper))
           }
           pvblurb <- paste("Monte Carlo test based on",
                            nsim, "simulations of CSR with fixed n")
         })

  statistic <- as.numeric(statistic)
  names(statistic) <- "A"
  
  out <- list(statistic=statistic,
              p.value=p.value,
              alternative=altblurb,
              method=c("Hopkins-Skellam test of CSR", pvblurb),
              data.name=Xname)
  class(out) <- "htest"
  return(out)
}
