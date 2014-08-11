clarkevans <- function(X, correction=c("none", "Donnelly", "cdf"),
                       clipregion=NULL)
{
  verifyclass(X, "ppp")
  W <- X$window

  # validate correction argument
  gavecorrection <- !missing(correction)
  correction <- pickoption("correction", correction,
                           c(none="none",
                             Donnelly="Donnelly",
                             donnelly="Donnelly",
                             guard="guard",
                             cdf="cdf"),
                           multi=TRUE)

  if(("Donnelly" %in% correction) && (W$type != "rectangle")) {
    if(gavecorrection)
      warning("Donnelly correction only available for rectangular windows")
    correction <- correction[correction != "Donnelly"]
  }

  # guard correction applied iff `clipregion' is present
  isguard <- "guard" %in% correction
  askguard <- any(isguard)
  gaveguard <- !is.null(clipregion)
  if(gaveguard)
    clipregion <- as.owin(clipregion)
  if(askguard && !gaveguard) {
    warning("guard correction not performed; clipregion not specified")
    correction <- correction[!isguard]
  } else if(gaveguard && !askguard) 
    correction <- c(correction, "guard")

  return(clarkevansCalc(X, correction, clipregion))
}

clarkevans.test <- function(X, ..., 
                            correction="none",
                            clipregion=NULL,
                            alternative=c("two.sided", "less", "greater"),
                            nsim=1000
                            ) {
  Xname <- short.deparse(substitute(X))
  verifyclass(X, "ppp")
  W <- X$window

  # validate SINGLE correction
  correction <- pickoption("correction", correction,
                           c(none="none",
                             Donnelly="Donnelly",
                             donnelly="Donnelly",
                             guard="guard",
                             cdf="cdf"))
  switch(correction,
         none={
           corrblurb <- "No edge correction"
         },
         Donnelly={
           if(W$type != "rectangle")
             stop("Donnelly correction only available for rectangular windows")
           corrblurb <- "Donnelly correction"
         },
         guard={
           if(is.null(clipregion))
             stop("clipregion not specified")
           clipregion <- as.owin(clipregion)
           corrblurb <- "Guard correction"
         },
         cdf={
           corrblurb <- "CDF correction"
         })

  # alternative hypothesis
  if(missing(alternative) || is.null(alternative))
    alternative <- "two.sided"
  alternative <- pickoption("alternative", alternative,
                           c(two.sided="two.sided",
                             less="less",
                             clustered="less",
                             greater="greater",
                             regular="greater"))

  altblurb <-
    switch(alternative,
           two.sided="two-sided",
           less="mean nn distance less than expected under CSR (clustered)",
           greater="mean nn distance greater than expected under CSR (regular)")

  # compute observed value
  statistic <- clarkevansCalc(X, correction=correction, clipregion=clipregion,
                              working=TRUE)
  working <- attr(statistic, "working")
  #
  if(correction == "none") {
    # standard Normal p-value
    SE <- with(working, sqrt(((4-pi)*area)/(4 * pi))/npts)
    Z <- with(working, (Dobs - Dpois)/SE)
    p.value <- switch(alternative,
                      less=pnorm(Z),
                      greater=1 - pnorm(Z),
                      two.sided= 2*(1-pnorm(abs(Z))))
    pvblurb <- "Z-test"
  } else {
    # Monte Carlo p-value
    sims <- numeric(nsim)
    intensity <- working$intensity
    for(i in 1:nsim) {
      Xsim <- rpoispp(intensity, win=W)
      sims[i] <- clarkevansCalc(Xsim, correction=correction,
                                clipregion=clipregion)
    }
    prob <- mean(sims <= statistic)
    p.value <- switch(alternative,
                      less=prob,
                      greater=1 - prob,
                      two.sided= 2*min(prob, 1-prob))
    
    pvblurb <- paste("Monte Carlo test based on",
                     nsim, "simulations of CSR")
  }

  statistic <- as.numeric(statistic)
  names(statistic) <- "R"
  
  out <- list(statistic=statistic,
              p.value=p.value,
              alternative=altblurb,
              method=c("Clark-Evans test", corrblurb, pvblurb),
              data.name=Xname)
  class(out) <- "htest"
  return(out)
}

clarkevansCalc <- function(X, correction="none", clipregion=NULL,
                           working=FALSE) {
  # calculations for Clark-Evans index or test
  W <- X$window
  area <- area.owin(W)
  npts <- npoints(X)
  intensity <- npts/area
  # R undefined for empty point pattern
  if(npts == 0)
    return(NA)
  # Dobs = observed mean nearest neighbour distance
  nndistX <- nndist(X)
  Dobs <- mean(nndistX)
  # Dpois = Expected mean nearest neighbour distance for Poisson process
  Dpois <- 1/(2*sqrt(intensity))

  statistic <- NULL
  
  # Naive uncorrected value
  if("none" %in% correction) {
    Rnaive <- Dobs/Dpois
    statistic <- c(statistic, naive=Rnaive)
  }
  # Donnelly edge correction
  if("Donnelly" %in% correction) {
     # Dedge = Edge corrected mean nearest neighbour distance, Donnelly 1978
    if(W$type == "rectangle") {
      perim <- perimeter(W)
      Dkevin  <- Dpois + (0.0514+0.0412/sqrt(npts))*perim/npts
      Rkevin <- Dobs/Dkevin
    } else 
      Rkevin <- NA
    statistic <- c(statistic, Donnelly=Rkevin)
  }
  # guard area method
  if("guard" %in% correction && !is.null(clipregion)) {
    # use nn distances from points inside `clipregion'
    ok <- inside.owin(X, , clipregion)
    Dguard <- mean(nndistX[ok])
    Rguard <- Dguard/Dpois
    statistic <- c(statistic, guard=Rguard)
  }
  if("cdf" %in% correction) {
    # compute mean of estimated nearest-neighbour distance distribution G
    G <- Gest(X)
    numer <- stieltjes(function(x){x}, G)$km
    denom <- stieltjes(function(x){rep(1, length(x))}, G)$km
    Dcdf <- numer/denom
    Rcdf <- Dcdf/Dpois
    statistic <- c(statistic, cdf=Rcdf)
  }
  if(working)
    attr(statistic, "working") <-
      list(area=area, npts=npts, intensity=intensity, Dpois=Dpois, Dobs=Dobs)

  return(statistic)
}
