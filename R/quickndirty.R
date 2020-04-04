#'
#'   quick-and-dirty KDE for points on a network
#'
#'   Copyright (C) 2019 Adrian Baddeley, Suman Rakshit and Tilman Davies
#'
#'   $Revision: 1.5 $ $Date: 2020/04/04 02:55:54 $

densityQuick.lpp <- function(X, sigma=NULL, ...,
                             kernel="gaussian",
                             at=c("pixels", "points"),
                             what=c("estimate", "se", "var"),
                             leaveoneout=TRUE,
                             diggle = FALSE,
                             edge2D = FALSE,
                             weights=NULL,
                             positive=FALSE) {
  #' kernel density estimation
  stopifnot(is.lpp(X))
  what <- match.arg(what)
  if(is.function(sigma)) sigma <- sigma(X)
  qkdeEngine(X=X, sigma=sigma, kernel=kernel,
             at=at, what=what, leaveoneout=leaveoneout,
             diggle=diggle, edge2D=edge2D,
             weights=weights, positive=positive,
             ...)
}

qkdeEngine <- function(X, sigma=NULL, ...,
                       at=c("pixels", "points"),
                       what=c("estimate", "se", "var"),
                       leaveoneout=TRUE,
                       diggle = FALSE,
                       raw=FALSE,
                       edge2D = FALSE,
                       edge = edge2D,
                       weights=NULL,
                       varcov=NULL,
                       positive=FALSE,
                       shortcut=TRUE,
                       precomputed=NULL,
                       savecomputed=FALSE) {
  stopifnot(is.lpp(X))
  at <- match.arg(at)
  what <- match.arg(what)
  L <- domain(X)
  S <- as.psp(L)
  XX <- as.ppp(X)
  stuff <- resolve.2D.kernel(x=XX, sigma=sigma, varcov=varcov, ...)
  sigma <- stuff$sigma
  varcov <- stuff$varcov

  if(is.infinite(stuff$cutoff)) {
    #' infinite bandwidth
    result <-
      switch(at,
             pixels = {
               as.linim(flatdensityfunlpp(X, weights=weights, what=what))
             },
             points = {
               flatdensityatpointslpp(X, weights=weights, what=what)
             })
    attr(result, "sigma") <- sigma
    attr(result, "varcov") <- varcov
    return(result)
  }
                     
  switch(what,
         estimate = {
           if(shortcut) {
             PS <- precomputed$PS %orifnull% pixellate(S, ...,
                                                       DivideByPixelArea=TRUE)
             KS <- blur(PS, sigma, normalise=edge2D, bleed=FALSE,
                        ..., varcov=varcov)
           } else {
             KS <- density(S, sigma, ..., edge=edge2D, varcov=varcov)
           }
           if(diggle && !raw) 
             weights <- (weights %orifnull% 1) / KS[XX]
           KX <- density(XX, sigma, ..., weights=weights,
                         at=at, leaveoneout=leaveoneout,
                         edge=edge2D, diggle=FALSE, positive=FALSE,
                         varcov=varcov)
         },
         se = ,
         var= {
           tau <- taumat <- NULL
           if(is.null(varcov)) {
             varconst <- 1/(4 * pi * prod(ensure2vector(sigma)))
             tau <- sigma/sqrt(2)
           } else {
             varconst <- 1/(4 * pi * sqrt(det(varcov)))
             taumat <- varcov/2
           }
           if(shortcut) {
             PS <- precomputed$PS %orifnull% pixellate(S, ...,
                                                       DivideByPixelArea=TRUE)
             KS <- blur(PS, sigma,
                        normalise=edge2D, bleed=FALSE, varcov=varcov)^2
           } else {
             KS <- density(S, sigma, ..., edge=edge2D, varcov=varcov)^2
           }
           if(diggle && !raw) 
             weights <- (weights %orifnull% 1) / KS[XX]
           KX <- varconst * density(XX, sigma=tau, ..., weights=weights,
                                 at=at, leaveoneout=leaveoneout,
                                 edge=edge2D, diggle=FALSE, positive=FALSE,
                                 varcov=taumat)
         })
  switch(at,
         points = {
           result <- if(diggle || raw) KX else (KX/(KS[XX]))
           if(positive)
             result <- pmax(result, .Machine$double.xmin)
           if(savecomputed) {
             #' save geometry info for re-use
             savedstuff <- list(PS=PS, M=solutionset(PS > 0), df=NULL)
           }
         },
         pixels = {
           Z <- if(diggle || raw) KX else (KX/KS)
           M <- if(shortcut) {
                  precomputed$M %orifnull% solutionset(PS > 0)
                } else as.mask.psp(S, KS)
           Z <- Z[M, drop=FALSE]
           #' build linim object, using precomputed sample points if available
           result <- linim(L, Z, restrict=FALSE, df=precomputed$df)
           if(positive)
             result <- eval.linim(pmax(result, .Machine$double.xmin))
           if(savecomputed) {
             #' save geometry info for re-use
             dfg <- attr(result, "df")
             dfg <- dfg[, colnames(dfg) != "values"]
             savedstuff <- list(PS=PS, M=M, df=dfg)
           }
         })
  if(what == "se") result <- sqrt(result)
  attr(result, "sigma") <- sigma
  attr(result, "varcov") <- varcov
  if(raw) attr(result, "denominator") <- KS
  if(savecomputed) attr(result, "savedstuff") <- savedstuff
  return(result)
}

