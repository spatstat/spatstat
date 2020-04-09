#
#    relrisk.lpp.R
#
#   Estimation of relative risk on network
#
#  $Revision: 1.2 $  $Date: 2020/04/09 03:46:14 $
#

relrisk.lpp <- local({

  relrisk.lpp <- function(X, sigma, ..., 
                          at=c("pixels", "points"), 
                          relative=FALSE,
                          adjust=1, 
                          casecontrol=TRUE, control=1, case,
                          finespacing=FALSE) {
    stopifnot(is.lpp(X))
    stopifnot(is.multitype(X))
    control.given <- !missing(control)
    case.given <- !missing(case)
    at <- match.arg(at)
    npts <- npoints(X)
    Y <- split(X)
    uX <- unmark(X)
    types <- names(Y)
    ntypes <- length(Y)
    if(ntypes == 1)
      stop("Data contains only one type of points")
    casecontrol <- casecontrol && (ntypes == 2)
    if((control.given || case.given) && !(casecontrol || relative)) {
      aa <- c("control", "case")[c(control.given, case.given)]
      nn <- length(aa)
      warning(paste(ngettext(nn, "Argument", "Arguments"),
                    paste(sQuote(aa), collapse=" and "),
                    ngettext(nn, "was", "were"),
                    "ignored, because relative=FALSE and",
                    if(ntypes==2) "casecontrol=FALSE" else
                    "there are more than 2 types of points"))
    }
    marx <- marks(X)
    imarks <- as.integer(marx)
    lev <- levels(marx)
    ## compute bandwidth (no bandwidth selection yet)
    check.1.real(sigma) # could be Inf
    sigma <- adjust * sigma
    ## .........................................
    ## compute intensity estimates for each type
    ## .........................................
    switch(at,
           pixels = {
             ## intensity estimates of each type
             Deach <- solapply(Y, density.lpp, sigma=sigma,
                               ..., finespacing=finespacing)
             ## compute intensity estimate for unmarked pattern
             Dall  <- density(X, sigma=sigma,
                              ..., finespacing=finespacing)
           },
           points = {
             ## intensity estimates of each type **at each data point**
             Deachfun <- solapply(Y, densityfun.lpp, sigma=sigma,
                                  ..., finespacing=finespacing)
             Deach <- as.data.frame(sapply(Deachfun, function(f, P) f(P), P=X))
             ## leave-one-out estimates
             Dself <- lapply(Y, density.lpp, sigma=sigma,
                             at="points", leaveoneout=TRUE,
                             ..., finespacing=finespacing)
             ## insert leave-one-out estimates in correct place
             Deachsplit <- split(Deach, marx)
             for(j in 1:ntypes) {
               Deachsplit[[j]][, j] <- Dself[[j]]
             }
             split(Deach, marx) <- Deachsplit
             ## total
             Dall <- rowSums(Deach)
           })
    ## .........................................
    ## compute probabilities/risks
    ## .........................................
    if(ntypes == 2 && casecontrol) {
      if(control.given || !case.given) {
        stopifnot(length(control) == 1)
        if(is.numeric(control)) {
          icontrol <- control <- as.integer(control)
          stopifnot(control %in% 1:2)
        } else if(is.character(control)) {
          icontrol <- match(control, levels(marks(X)))
          if(is.na(icontrol)) stop(paste("No points have mark =", control))
        } else
          stop(paste("Unrecognised format for argument", sQuote("control")))
        if(!case.given)
          icase <- 3 - icontrol
      }
      if(case.given) {
        stopifnot(length(case) == 1)
        if(is.numeric(case)) {
          icase <- case <- as.integer(case)
          stopifnot(case %in% 1:2)
        } else if(is.character(case)) {
          icase <- match(case, levels(marks(X)))
          if(is.na(icase)) stop(paste("No points have mark =", case))
        } else stop(paste("Unrecognised format for argument", sQuote("case")))
        if(!control.given) 
          icontrol <- 3 - icase
      }
      ## compute ......
      switch(at,
             pixels = {
               ## compute probability of case
               pcase <- Deach[[icase]]/Dall
               ## correct small numerical errors
               pcase <- clamp01(pcase)
               ## trap NaN values
               nbg <- badvalues(pcase)
               if(any(nbg)) {
                 ## apply l'Hopital's rule:
                 ##     p(case) = 1{nearest neighbour is case}
                 distcase <- as.linim(distfun(Y[[icase]]))
                 distcontrol <- as.linim(distfun(Y[[icontrol]]))
                 closecase <- eval.linim(as.integer(distcase < distcontrol))
                 pcase[nbg] <- closecase[nbg]
               }
               if(!relative) {
                 result <- pcase
               } else {
                 result <- eval.im(ifelse(pcase < 1, pcase/(1-pcase), NA))
               }
             },
             points={
               ## compute probability of case
               pcase <- Deach[,icase]/Dall
               ## correct small numerical errors
               pcase <- clamp01(pcase)
               ## trap NaN values
               if(any(nbg <- badvalues(pcase))) {
                 ## apply l'Hopital's rule
                 nntype <- imarks[nnwhich(X)]
                 pcase[nbg] <- as.integer(nntype[nbg] == icase)
               }
               if(!relative) {
                 result <- pcase
               } else {
                 result <- ifelse(pcase < 1, pcase/(1-pcase), NA)
               }
             })
    } else {
      ## several types
      if(relative) {
        ## need 'control' type
        stopifnot(length(control) == 1)
        if(is.numeric(control)) {
          icontrol <- control <- as.integer(control)
          stopifnot(control %in% 1:ntypes)
        } else if(is.character(control)) {
          icontrol <- match(control, levels(marks(X)))
          if(is.na(icontrol)) stop(paste("No points have mark =", control))
        } else
          stop(paste("Unrecognised format for argument", sQuote("control")))
      }
      switch(at,
             pixels={
               probs <- as.solist(lapply(Deach, "/", e2=Dall))
               ## correct small numerical errors
               probs <- as.solist(lapply(probs, clamp01))
               ## trap NaN values
               nbg <- lapply(probs, badvalues)
               nbg <- Reduce("|", nbg)
               if(any(nbg)) {
                 ## apply l'Hopital's rule
                 distX <- as.linim(distfun(X))
                 whichnn <- as.linim(nnfun(X))
                 typenn <- eval.im(imarks[whichnn])
                 typennsub <- typenn[nbg]
                 for(k in seq_along(probs)) 
                   probs[[k]][nbg] <- (typennsub == k)
               }
               if(!relative) {
                 result <- probs
               } else {
                 result <- solapply(probs,
                                    divideifpositive,
                                    d = probs[[icontrol]])
               }
             },
             points = {
               probs <- Deach/Dall
               ## correct small numerical errors
               probs <- clamp01(probs)
               ## trap NaN values
               bad <- badvalues(probs)
               badrow <- matrowany(bad)
               if(any(badrow)) {
                 ## apply l'Hopital's rule
                 typenn <- imarks[nnwhich(X)]
                 probs[badrow, ] <- (typenn == col(result))[badrow, ]
               }
               if(!relative) {
                 result <- probs
               } else {
                 result <- probs/probs[,icontrol]
               }
            })
    }
    attr(result, "sigma") <- sigma
    return(result)
  }

  clamp01 <- function(x) {
    if(is.linim(x)) return(eval.linim(pmin(pmax(x, 0), 1)))
    if(is.im(x)) return(eval.im(pmin(pmax(x, 0), 1)))
    if(is.data.frame(x)) x <- as.matrix(x)
    return(pmin(pmax(x, 0), 1))
  }

  badvalues <- function(x) {
    if(is.linim(x)) return(eval.linim(!is.finite(x)))
    if(is.im(x)) return(eval.im(!is.finite(x)))
    if(is.data.frame(x)) x <- as.matrix(x)
    return(!(is.finite(x) | is.na(x)))
  }

  divideifpositive <- function(z, d) { eval.linim(ifelse(d > 0, z/d, NA)) }
  
  relrisk.lpp
})


