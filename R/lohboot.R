#
#  lohboot.R
#
#  $Revision: 1.7 $   $Date: 2013/08/01 09:40:03 $
#
#  Loh's bootstrap CI's for local pcf, local K etc
#

lohboot <-
  function(X,
           fun=c("pcf", "Kest", "Lest", "pcfinhom", "Kinhom", "Linhom"),
           ..., nsim=200, confidence=0.95, global=FALSE, type=7) {
  stopifnot(is.ppp(X))
  fun <- match.arg(fun)
  # validate confidence level
  stopifnot(confidence > 0.5 && confidence < 1)
  alpha <- 1 - confidence
  if(!global) {
    probs <- c(alpha/2, 1-alpha/2)
    rank <- nsim * probs[2]
  } else {
    probs <- 1-alpha
    rank <- nsim * probs
  }
  if(abs(rank - round(rank)) > 0.001)
    warning(paste("confidence level", confidence,
                  "corresponds to a non-integer rank", paren(rank),
                  "so quantiles will be interpolated"))
  n <- npoints(X)
  # compute local functions
  localfun <- switch(fun,
                     pcf=localpcf,
                     Kest=localK,
                     Lest=localK,
                     pcfinhom=localpcfinhom,
                     Kinhom=localKinhom,
                     Linhom=localKinhom)
  f <- localfun(X, ...)
  theo <- f$theo
  # parse edge correction info
  correction <- attr(f, "correction")
  switch(correction,
         none      = { ctag <- "un";    cadj <- "uncorrected" },
         border    = { ctag <- "bord";  cadj <- "border-corrected" },
         translate = { ctag <- "trans"; cadj <- "translation-corrected" },
         isotropic = { ctag <- "iso";   cadj <- "Ripley isotropic corrected" })
  # first n columns are the local pcfs (etc) for the n points of X
  y <- as.matrix(as.data.frame(f))[, 1:n]
  # average them
  ymean <- rowMeans(y, na.rm=TRUE)
  # resample
  ystar <- matrix(, nrow=nrow(y), ncol=nsim)
  for(i in 1:nsim) {
    # resample n points with replacement
    ind <- sample(n, replace=TRUE)
    # average their local pcfs
    ystar[,i] <- rowMeans(y[,ind], na.rm=TRUE)
  }
  # compute quantiles
  if(!global) {
    # pointwise quantiles
    hilo <- apply(ystar, 1, quantile,
                  probs=probs, na.rm=TRUE, type=type)
  } else {
    # quantiles of deviation
    ydif <- sweep(ystar, 1, ymean)
    ydev <- apply(abs(ydif), 2, max, na.rm=TRUE)
    crit <- quantile(ydev, probs=probs, na.rm=TRUE, type=type)
    hilo <- rbind(ymean - crit, ymean + crit)
  }
  # now transform from K to L if required
  if(fun %in% c("Lest", "Linhom")) {
    ymean <- sqrt(ymean/pi)
    theo  <- sqrt(theo/pi)
    hilo  <- sqrt(hilo/pi)
  }
  # create fv object
  df <- data.frame(r=f$r,
                   theo=theo,
                   ymean,
                   lo=hilo[1,],
                   hi=hilo[2,])
  colnames(df)[3] <- ctag
  CIlevel <- paste(100 * confidence, "%% confidence", sep="")
  desc <- c("distance argument r",
            "theoretical Poisson %s",
            paste(cadj, "estimate of %s"),
            paste("lower", CIlevel, "limit for %s"),
            paste("upper", CIlevel, "limit for %s"))
  clabl <- paste("hat(%s)[", ctag, "](r)", sep="")
  labl <- c("r", "%s[pois](r)", clabl, "%s[loCI](r)", "%s[hiCI](r)")
  switch(fun,
         pcf={ fname <- "g" ; ylab <- quote(g(r)) },
         Kest={ fname <- "K" ; ylab <- quote(K(r)) },
         Lest={ fname <- "L" ; ylab <- quote(L(r)) },
         pcfinhom={ fname <- "g[inhom]" ; ylab <- quote(g[inhom](r)) },
         Kinhom={ fname <- "K[inhom]" ; ylab <- quote(K[inhom](r)) },
         Linhom={ fname <- "L[inhom]" ; ylab <- quote(L[inhom](r)) })
  g <- fv(df, "r", ylab, ctag, , c(0, max(f$r)), labl, desc, fname=fname)
  formula(g) <- . ~ r
  fvnames(g, ".") <- c(ctag, "hi", "lo", "theo")
  fvnames(g, ".s") <- c("hi", "lo")
  unitname(g) <- unitname(X)
  g
}


    
  
  
  
  
