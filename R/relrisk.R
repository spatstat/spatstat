#
#    relrisk.R
#
#   Estimation of relative risk
#
#  $Revision: 1.31 $  $Date: 2016/09/25 10:22:46 $
#

relrisk <- function(X, ...) UseMethod("relrisk")
                                      
relrisk.ppp <- local({

  relrisk.ppp <- function(X, sigma=NULL, ..., varcov=NULL, at="pixels",
                      relative=FALSE, se=FALSE,
                      casecontrol=TRUE, control=1, case) {
    stopifnot(is.ppp(X))
    stopifnot(is.multitype(X))
    control.given <- !missing(control)
    case.given <- !missing(case)
    if(!relative && (control.given || case.given)) {
      aa <- c("control", "case")[c(control.given, case.given)]
      nn <- length(aa)
      warning(paste(ngettext(nn, "Argument", "Arguments"),
                    paste(sQuote(aa), collapse=" and "),
                    ngettext(nn, "was", "were"),
                    "ignored, because relative=FALSE"))
    }
    npts <- npoints(X)
    Y <- split(X)
    uX <- unmark(X)
    types <- names(Y)
    ntypes <- length(Y)
    if(ntypes == 1)
      stop("Data contains only one type of points")
    marx <- marks(X)
    imarks <- as.integer(marx)
    lev <- levels(marx)
    ## trap arguments
    dotargs <- list(...)
    isbwarg <- names(dotargs) %in% c("method", "nh", "hmin", "hmax", "warn")
    bwargs <- dotargs[isbwarg]
    dargs  <- dotargs[!isbwarg]
    ## using edge corrections?
    edge   <- resolve.1.default(list(edge=TRUE), list(...))
    diggle <- resolve.1.default(list(diggle=FALSE), list(...))
    ## bandwidth
    if(is.null(sigma) && is.null(varcov)) 
      sigma <- do.call(bw.relrisk, append(list(X), bwargs))
    SmoothPars <- append(list(sigma=sigma, varcov=varcov, at=at), dargs)
    if(se) {
      ## determine other bandwidth for variance estimation
      if(is.null(varcov)) {
        varconst <- 1/(4 * pi * prod(sigma))
        VarPars <- append(list(sigma=sigma/sqrt(2), at=at), dargs)
      } else {
        varconst <- 1/(4 * pi * sqrt(det(varcov)))
        VarPars <- append(list(varcov=varcov/2, at=at), dargs)
      }
      if(edge) {
        ## evaluate edge correction weights
        edgeim <- second.moment.calc(uX, sigma, what="edge", ...,
                                     varcov=varcov)
        if(diggle || at == "points") {
          edgeX <- safelookup(edgeim, uX, warn=FALSE)
          diggleX <- 1/edgeX
          diggleX[!is.finite(diggleX)] <- 0
        }
        edgeim <- edgeim[Window(X), drop=FALSE]
      }
    }
    ## .........................................
    ## compute intensity estimates for each type
    ## .........................................
    switch(at,
           pixels = {
             ## intensity estimates of each type
             Deach <- do.call(density.splitppp,
                              append(list(x=Y), SmoothPars))
             ## compute intensity estimate for unmarked pattern
             Dall <- Reduce("+", Deach)
             ## variance terms
             if(se) {
               if(!edge) {
                 ## no edge correction
                 Veach <- do.call(density.splitppp,
                                  append(list(x=Y), VarPars))
               } else if(!diggle) {
                 ## edge correction e(u)
                 Veach <- do.call(density.splitppp,
                                  append(list(x=Y), VarPars))
                 Veach <- lapply(Veach, "/", e2=edgeim)
               } else {
                 ## Diggle edge correction e(x_i)
                 Veach <- mapply(density.ppp,
                                 x=Y,
                                 weights=split(diggleX, marx),
                                 MoreArgs=VarPars,
                                 SIMPLIFY=FALSE)
               }
               Veach <- lapply(Veach, "*", varconst)
               Vall <- Reduce("+", Veach)
             }
           },
           points = {
             ## intensity estimates of each type **at each data point**
             ## dummy variable matrix
             dumm <- matrix(0, npts, ntypes)
             dumm[cbind(seq_len(npts), imarks)] <- 1
             colnames(dumm) <- lev
             Deach <- do.call(density.ppp,
                              append(list(x=uX, weights=dumm),
                                     SmoothPars))
             ## compute intensity estimate for unmarked pattern
             Dall <- rowSums(Deach)
             ## variance terms
             if(se) {
               if(!edge) {
                 ## no edge correction
                 Veach <- do.call(density.ppp,
                                  append(list(x=uX, weights=dumm),
                                         VarPars))
               } else if(!diggle) {
                 ## edge correction e(u)
                 Veach <- do.call(density.ppp,
                                  append(list(x=uX, weights=dumm),
                                         VarPars))
                 Veach <- Veach * diggleX
               } else {
                 ## Diggle edge correction e(x_i)
                 Veach <- do.call(density.ppp,
                                  append(list(x=uX, weights=dumm * diggleX),
                                         VarPars))
               }
               Veach <- Veach * varconst
               Vall <- rowSums(Veach)
             }
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
                 distcase <- distmap(Y[[icase]], xy=pcase)
                 distcontrol <- distmap(Y[[icontrol]], xy=pcase)
                 closecase <- eval.im(as.integer(distcase < distcontrol))
                 pcase[nbg] <- closecase[nbg]
               }
               if(!relative) {
                 if(!se) {
                   result <- pcase
                 } else {
                   Vcase <- Veach[[icase]]
                   NUM <- eval.im(Vcase * (1-2*pcase) + Vall * pcase^2)
                   SE <- eval.im(sqrt(pmax(NUM, 0))/Dall)
                   result <- list(estimate=pcase, SE=SE)
                 }
               } else {
                 rcase <- eval.im(ifelse(pcase < 1, pcase/(1-pcase), NA))
                 if(!se) {
                   result <- rcase
                 } else {
                   Vcase <- Veach[[icase]]
                   Vctrl <- Veach[[icontrol]]
                   Dctrl <- Deach[[icontrol]]
                   NUM <- eval.im(Vcase + Vctrl * rcase^2)
                   SE <- eval.im(sqrt(pmax(NUM, 0))/Dctrl)
                   result <- list(estimate=rcase, SE=SE)
                 }
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
                 if(!se) {
                   result <- pcase
                 } else {
                   NUM <- Veach[,icase] * (1-2*pcase) + Vall * pcase^2
                   SE <- sqrt(pmax(NUM, 0))/Dall
                   result <- list(estimate=pcase, SE=SE)
                 }
               } else {
                 rcase <- ifelse(pcase < 1, pcase/(1-pcase), NA)
                 if(!se) {
                   result <- rcase
                 } else {
                   NUM <- Veach[,icase] + Veach[,icontrol] * rcase^2
                   SE <- sqrt(pmax(NUM, 0))/Deach[,icontrol]
                   result <- list(estimate=rcase, SE=SE)
                 }
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
                 distX <- distmap(X, xy=Dall)
                 whichnn <- attr(distX, "index")
                 typenn <- eval.im(imarks[whichnn])
                 typennsub <- as.matrix(typenn)[nbg]
                 for(k in seq_along(result)) 
                   probs[[k]][nbg] <- (typennsub == k)
               }
               if(!relative) {
                 if(!se) {
                   result <- probs
                 } else {
                   SE <- list()
                   for(i in 1:ntypes) {
                     NUM <- (Veach[[i]] * (1 - 2 * probs[[i]])
                             + Vall * probs[[i]]^2)
                     SE[[i]] <- eval.im(sqrt(pmax(NUM, 0))/Dall)
                   }
                   SE <- as.solist(SE)
                   names(SE) <- types
                   result <- list(estimate=probs, SE=SE)
                 }
               } else {
                 risks <- as.solist(lapply(probs,
                                           function(z, d) {
                                             eval.im(ifelse(d > 0, z/d, NA))
                                           },
                                           d = probs[[icontrol]]))
                 if(!se) {
                   result <- risks
                 } else {
                   Vctrl <- Veach[[icontrol]]
                   Dctrl <- Deach[[icontrol]]
                   SE <- list()
                   for(i in 1:ntypes) {
                     NUM <- Veach[[i]] + Vctrl * risks[[i]]^2
                     SE[[i]] <- eval.im(sqrt(pmax(NUM, 0))/Dctrl)
                   }
                   SE <- as.solist(SE)
                   names(SE) <- types
                   result <- list(estimate=risks, SE=SE)
                   
                 }
               }
             },
             points = {
               probs <- Deach/Dall
               ## correct small numerical errors
               probs <- clamp01(probs)
               ## trap NaN values
               bad <- badvalues(probs)
               badrow <- apply(bad, 1, any)
               if(any(badrow)) {
                 ## apply l'Hopital's rule
                 typenn <- imarks[nnwhich(X)]
                 probs[badrow, ] <- (typenn == col(result))[badrow, ]
               }
               if(!relative) {
                 if(!se) {
                   result <- probs
                 } else {
                   NUM <- Veach * (1-2*probs) + Vall * probs^2
                   SE <- sqrt(pmax(NUM, 0))/Dall
                   result <- list(estimate=probs, SE=SE)
                }
               } else {
                 risks <- probs/probs[,icontrol]
                 if(!se) {
                   result <- risks
                 } else {
                   NUM <- Veach + Veach[,icontrol] * risks^2
                   NUM[,icontrol] <- 0
                   SE <- sqrt(pmax(NUM, 0))/Deach[,icontrol]
                   result <- list(estimate=risks, SE=SE)
                 }
               }
            })
    }
    attr(result, "sigma") <- sigma
    attr(result, "varcov") <- varcov
    return(result)
  }

  clamp01 <- function(x) {
    if(is.im(x)) return(eval.im(pmin(pmax(x, 0), 1)))
    return(pmin(pmax(x, 0), 1))
  }

  badvalues <- function(x) {
    if(is.im(x)) x <- as.matrix(x)
    return(!(is.finite(x) | is.na(x)))
  }

  reciprocal <- function(x) 1/x
  
  relrisk.ppp
})


bw.stoyan <- function(X, co=0.15) {
  ## Stoyan's rule of thumb
  stopifnot(is.ppp(X))
  n <- npoints(X)
  W <- Window(X)
  a <- area(W)
  stoyan <- co/sqrt(5 * n/a)
  return(stoyan)
}


bw.relrisk <- function(X, method="likelihood",
                       nh=spatstat.options("n.bandwidth"),
                       hmin=NULL, hmax=NULL, warn=TRUE) {
  stopifnot(is.ppp(X))
  stopifnot(is.multitype(X))
  ## rearrange in ascending order of x-coordinate (for C code)
  X <- X[fave.order(X$x)]
  ##
  Y <- split(X)
  ntypes <- length(Y)
  if(ntypes == 1)
    stop("Data contains only one type of points")
  marx <- marks(X)
  method <- pickoption("method", method,
                       c(likelihood="likelihood",
                         leastsquares="leastsquares",
                         ls="leastsquares",
                         LS="leastsquares",
                         weightedleastsquares="weightedleastsquares",
                         wls="weightedleastsquares",
                         WLS="weightedleastsquares"))
  ## 
  if(method != "likelihood") {
    ## dummy variables for each type
    imarks <- as.integer(marx)
    if(ntypes == 2) {
      ## 1 = control, 2 = case
      indic <- (imarks == 2)
      y01   <- as.integer(indic)
    } else {
      indic <- matrix(FALSE, n, ntypes)
      indic[cbind(seq_len(n), imarks)] <- TRUE
      y01  <- indic * 1
    }
    X01 <- X %mark% y01
  }
  ## cross-validated bandwidth selection
  ## determine a range of bandwidth values
  n <- npoints(X)
  if(is.null(hmin) || is.null(hmax)) {
    W <- Window(X)
    a <- area(W)
    d <- diameter(as.rectangle(W))
    ## Stoyan's rule of thumb applied to the least and most common types
    mcount <- table(marx)
    nmin <- max(1, min(mcount))
    nmax <- max(1, max(mcount))
    stoyan.low <- 0.15/sqrt(nmax/a)
    stoyan.high <- 0.15/sqrt(nmin/a)
    if(is.null(hmin)) 
      hmin <- max(minnndist(unique(X)), stoyan.low/5)
    if(is.null(hmax)) {
      hmax <- min(d/4, stoyan.high * 20)
      hmax <- max(hmax, hmin * 2)
    }
  } else stopifnot(hmin < hmax)
  ##
  h <- geomseq(from=hmin, to=hmax, length.out=nh)
  cv <- numeric(nh)
  ## 
  ## compute cross-validation criterion
  switch(method,
         likelihood={
           methodname <- "Likelihood"
           ## for efficiency, only compute the estimate of p_j(x_i)
           ## when j = m_i = mark of x_i.
           Dthis <- numeric(n)
           for(i in seq_len(nh)) {
             Dall <- density.ppp(X, sigma=h[i], at="points", edge=FALSE,
                                 sorted=TRUE)
             Deach <- density.splitppp(Y, sigma=h[i], at="points", edge=FALSE,
                                       sorted=TRUE)
             split(Dthis, marx) <- Deach
             pthis <- Dthis/Dall
             cv[i] <- -mean(log(pthis))
           }
         },
         leastsquares={
           methodname <- "Least Squares"
           for(i in seq_len(nh)) {
             phat <- Smooth(X01, sigma=h[i], at="points", leaveoneout=TRUE,
                            sorted=TRUE)
             cv[i] <- mean((y01 - phat)^2)
           }
         },
         weightedleastsquares={
           methodname <- "Weighted Least Squares"
           ## need initial value of h from least squares
           h0 <- bw.relrisk(X, "leastsquares", nh=ceiling(nh/4))
           phat0 <- Smooth(X01, sigma=h0, at="points", leaveoneout=TRUE,
                           sorted=TRUE)
           var0 <- phat0 * (1-phat0)
           var0 <- pmax.int(var0, 1e-6)
           for(i in seq_len(nh)) {
             phat <- Smooth(X01, sigma=h[i], at="points", leaveoneout=TRUE,
                            sorted=TRUE)
             cv[i] <- mean((y01 - phat)^2/var0)
           }
         })
  ## optimize
  iopt <- which.min(cv)
  ##
  if(warn && (iopt == nh || iopt == 1)) 
    warning(paste("Cross-validation criterion was minimised at",
                  if(iopt == 1) "left-hand" else "right-hand",
                  "end of interval",
                  "[", signif(hmin, 3), ",", signif(hmax, 3), "];",
                  "use arguments hmin, hmax to specify a wider interval"))
  ##    
  result <- bw.optim(cv, h, iopt,
                     hname="sigma", 
                     creator="bw.relrisk",
                     criterion=paste(methodname, "Cross-Validation"))
  return(result)
}

which.max.im <- function(x) {
  .Deprecated("im.apply", "spatstat",
              "which.max.im(x) is deprecated: use im.apply(x, which.max)")
  ans <- im.apply(x, which.max)
  return(ans)
}

