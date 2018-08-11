#
#  ho.R
#
#  Huang-Ogata method 
#
#  $Revision: 1.17 $ $Date: 2016/03/15 07:42:26 $
#

ho.engine <- function(model, ..., nsim=100, nrmh=1e5,
                        start=NULL,
                        control=list(nrep=nrmh), verb=TRUE) {
  verifyclass(model, "ppm")

  if(is.null(start)) 
    start <- list(n.start=data.ppm(model)$n)
  
  # check that the model can be simulated
  if(!valid.ppm(model)) {
    warning("Fitted model is invalid - cannot be simulated")
    return(NULL)
  }
  
  # compute the observed value of the sufficient statistic
  X <- data.ppm(model)
  sobs <- suffstat(model, X)
  
  # generate 'nsim' realisations of the fitted model
  # and compute the sufficient statistics of the model
  rmhinfolist <- rmh(model, start, control, preponly=TRUE, verbose=FALSE)
  if(verb) {
    cat("Simulating... ")
    state <- list()
  }
  ndone <- 0
  while(ndone < nsim) {
    Xi <- rmhEngine(rmhinfolist, verbose=FALSE)
    v <- try(suffstat(model,Xi))
    if(!inherits(v, "try-error")) {
      if(ndone == 0) 
        svalues <- matrix(, nrow=nsim, ncol=length(v))
      ndone <- ndone + 1
      svalues[ndone, ] <- v
    }
    if(verb) state <- progressreport(ndone, nsim, state=state)
  }
  if(verb) cat("Done.\n\n")
  # calculate the sample mean and variance of the
  # sufficient statistic for the simulations
  smean <- apply(svalues, 2, mean, na.rm=TRUE)
  svar <- var(svalues, na.rm=TRUE)
  # value of canonical parameter from MPL fit
  theta0 <- coef(model)
  # Newton-Raphson update
  Vinverse <- solve(svar)
  theta <- theta0 + as.vector(Vinverse %*% (sobs - smean))
  ## appropriate names
  nama <- names(theta0)
  if(!is.null(nama)) {
    names(theta) <- nama
    dimnames(svar) <- dimnames(Vinverse) <- list(nama, nama)
  }
  ## update model
  newmodel <- model
  newmodel$coef <- theta
  newmodel$coef.orig <- theta0
  newmodel$method <- "ho"
  newmodel$fitter <- "ho"
  newmodel$fisher <- svar
  newmodel$varcov <- Vinverse
  # recompute fitted interaction
  newmodel$fitin <- NULL
  newmodel$fitin <- fitin(newmodel)
  ## update pseudolikelihood value using code in logLik.ppm
  newmodel$maxlogpl.orig <- model$maxlogpl
  newmodel$maxlogpl <- logLik(newmodel, new.coef=theta, warn=FALSE)
  ##
  return(newmodel)
}

