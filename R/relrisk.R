#
#    relrisk.R
#
#   Estimation of relative risk
#
#  $Revision: 1.17 $  $Date: 2013/02/17 23:48:03 $
#

relrisk <- function(X, sigma=NULL, ..., varcov=NULL, at="pixels",
                    casecontrol=TRUE) {
  stopifnot(is.ppp(X))
  stopifnot(is.multitype(X))
  Y <- split(X)
  ntypes <- length(Y)
  if(ntypes == 1)
    stop("Data contains only one type of points")
  marx <- marks(X)
  imarks <- as.integer(marx)
  lev <- levels(marx)
  # trap arguments
  dotargs <- list(...)
  isbwarg <- names(dotargs) %in% c("method", "nh", "hmin", "hmax", "warn")
  bwargs <- dotargs[isbwarg]
  dargs  <- dotargs[!isbwarg]
  # bandwidth
  if(is.null(sigma) && is.null(varcov)) {
    sigma <- do.call(bw.relrisk, append(list(X), bwargs))
  }
  # compute probabilities
  if(ntypes == 2 && casecontrol) {
    # 1 = control, 2 = case
    # compute densities
    Deach <- do.call(density.splitppp,
                     append(list(Y, sigma=sigma, varcov=varcov, at=at),
                            dargs))
    Dall <- do.call(density.ppp,
                    append(list(X, sigma=sigma, varcov=varcov, at=at),
                           dargs))
    # compute probability of case
    switch(at,
           pixels= {
             Dcase <- Deach[[2]]
             result <- eval.im(Dcase/Dall)
             # trap NaN values
             nbg <- as.matrix(eval.im(badprobability(result, FALSE)))
             if(any(nbg)) {
               # apply l'Hopital's rule:
               #     p(case) = 1{nearest neighbour is case}
               dist1 <- distmap(Y[[1]], xy=result)
               dist2 <- distmap(Y[[2]], xy=result)
               close2 <- eval.im(as.integer(dist2 < dist1))
               result[nbg] <- close2[nbg]
             }
           },
           points={
             result <- numeric(npoints(X))
             iscase <- (imarks == 2)
             result[iscase] <- Deach[[2]]/Dall[iscase]
             result[!iscase] <- 1 - Deach[[1]]/Dall[!iscase]
             # trap NaN values
             if(any(nbg <- badprobability(result, TRUE))) {
               # apply l'Hopital's rule
               nntype <- imarks[nnwhich(X)]
               result[nbg] <- as.integer(nntype[nbg] == 2)
             }
           })
  } else {
    # several types
    switch(at,
           pixels={
             Deach <- do.call(density.splitppp,
                              append(list(Y, sigma=sigma, varcov=varcov, at=at),
                                     dargs))
             Dall <- do.call(density.ppp,
                             append(list(X, sigma=sigma, varcov=varcov, at=at),
                                    dargs))
             result <- as.listof(lapply(Deach,
                                        function(d, dall) { eval.im(d/dall) },
                                        dall = Dall))
             # trap NaN values
             nbg <- lapply(result,
                           function(z) {
                             as.matrix(eval.im(badprobability(z, FALSE)))
                           })
             nbg <- Reduce("|", nbg)
             if(any(nbg)) {
               # apply l'Hopital's rule
               distX <- distmap(X, xy=Dall)
               whichnn <- attr(distX, "index")
               typenn <- eval.im(imarks[whichnn])
               typennsub <- as.matrix(typenn)[nbg]
               for(k in seq_along(result)) 
                 result[[k]][nbg] <- (typennsub == k)
             }
           },
           points = {
             npts <- npoints(X)
             # dummy variable matrix
             dumm <- matrix(0, npts, ntypes)
             dumm[cbind(seq_len(npts), imarks)] <- 1
             colnames(dumm) <- lev
             # compute probability of each type
             Z <- X %mark% dumm
             result <- do.call(smooth.ppp,
                               append(list(Z, sigma=sigma, varcov=varcov,
                                           at="points"),
                                      dargs))
             # trap NaN values
             bad <- badprobability(as.matrix(result), TRUE)
             badrow <- apply(bad, 1, any)
             if(any(badrow)) {
               # apply l'Hopital's rule
               typenn <- imarks[nnwhich(X)]
               result[badrow, ] <- (typenn == col(result))[badrow, ]
             }
           })
  }
  attr(result, "sigma") <- sigma
  attr(result, "varcov") <- varcov
  return(result)
}

bw.stoyan <- function(X, co=0.15) {
  # Stoyan's rule of thumb
  stopifnot(is.ppp(X))
  n <- npoints(X)
  W <- as.owin(X)
  a <- area.owin(W)
  stoyan <- co/sqrt(5 * n/a)
  return(stoyan)
}


bw.relrisk <- function(X, method="likelihood",
                       nh=spatstat.options("n.bandwidth"),
                       hmin=NULL, hmax=NULL, warn=TRUE) {
  stopifnot(is.ppp(X))
  stopifnot(is.multitype(X))
  # rearrange in ascending order of x-coordinate (for C code)
  X <- X[fave.order(X$x)]
  #
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
  # 
  if(method != "likelihood") {
    # dummy variables for each type
    imarks <- as.integer(marx)
    if(ntypes == 2) {
      # 1 = control, 2 = case
      indic <- (imarks == 2)
      y01   <- as.integer(indic)
    } else {
      indic <- matrix(FALSE, n, ntypes)
      indic[cbind(seq_len(n), imarks)] <- TRUE
      y01  <- indic * 1
    }
    X01 <- X %mark% y01
  }
  # cross-validated bandwidth selection
  # determine a range of bandwidth values
  n <- npoints(X)
  if(is.null(hmin) || is.null(hmax)) {
    W <- as.owin(X)
    a <- area.owin(W)
    d <- diameter(as.rectangle(W))
    # Stoyan's rule of thumb applied to the least and most common types
    mcount <- table(marx)
    nmin <- max(1, min(mcount))
    nmax <- max(1, max(mcount))
    stoyan.low <- 0.15/sqrt(nmax/a)
    stoyan.high <- 0.15/sqrt(nmin/a)
    if(is.null(hmin)) 
      hmin <- max(min(nndist(unique(X))), stoyan.low/5)
    if(is.null(hmax)) {
      hmax <- min(d/4, stoyan.high * 20)
      hmax <- max(hmax, hmin * 2)
    }
  } else stopifnot(hmin < hmax)
  #
  h <- exp(seq(from=log(hmin), to=log(hmax), length.out=nh))
  cv <- numeric(nh)
  # 
  # compute cross-validation criterion
  switch(method,
         likelihood={
           # for efficiency, only compute the estimate of p_j(x_i)
           # when j = m_i = mark of x_i.
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
           for(i in seq_len(nh)) {
             phat <- smooth.ppp(X01, sigma=h[i], at="points", leaveoneout=TRUE,
                                sorted=TRUE)
             cv[i] <- mean((y01 - phat)^2)
           }
         },
         weightedleastsquares={
           # need initial value of h from least squares
           h0 <- bw.relrisk(X, "leastsquares", nh=ceiling(nh/4))
           phat0 <- smooth.ppp(X01, sigma=h0, at="points", leaveoneout=TRUE,
                               sorted=TRUE)
           var0 <- phat0 * (1-phat0)
           var0 <- pmax(var0, 1e-6)
           for(i in seq_len(nh)) {
             phat <- smooth.ppp(X01, sigma=h[i], at="points", leaveoneout=TRUE,
                                sorted=TRUE)
             cv[i] <- mean((y01 - phat)^2/var0)
           }
         })
  # optimize
  iopt <- which.min(cv)
  #
  if(warn && (iopt == nh || iopt == 1)) 
    warning(paste("Cross-validation criterion was minimised at",
                  if(iopt == 1) "left-hand" else "right-hand",
                  "end of interval",
                  "[", signif(hmin, 3), ",", signif(hmax, 3), "];",
                  "use arguments hmin, hmax to specify a wider interval"))
  #    
  result <- bw.optim(cv, h, iopt,
                     xlab="sigma", ylab=paste(method, "CV"),
                     creator="bw.relrisk")
  return(result)
}

which.max.im <- function(x) {
  stopifnot(is.list(x))
  n <- length(x)
  if(n == 0)
    return(list())
  if(!all(unlist(lapply(x, is.im))))
    stop("x should be a list of images")
  nama <- names(x)
  xmax <- x[[1]]
  wmax <- eval.im(as.integer(xmax == xmax))
  if(n > 1) {
    for(i in 2:n) {
      xi <- x[[i]]
      xmaxnew <- eval.im(pmax(xi, xmax))
      wmaxnew <- eval.im(ifelse(xi > xmax, i, wmax))
      xmax <- xmaxnew
      wmax <- wmaxnew
    }
  }
  wmax <- eval.im(factor(wmax, levels=1:n))
  if(!is.null(nama))
    levels(wmax) <- nama
  return(wmax)
}
