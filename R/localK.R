#
#	localK.R		Getis-Franklin neighbourhood density function
#
#	$Revision: 1.21 $	$Date: 2015/07/11 08:19:26 $
#
#

"localL" <-
  function(X, ..., correction="Ripley", verbose=TRUE, rvalue=NULL)
{
  localK(X, wantL=TRUE,
         correction=correction, verbose=verbose, rvalue=rvalue)
}

"localLinhom" <-
  function(X, lambda=NULL, ..., correction="Ripley", verbose=TRUE, rvalue=NULL,
           sigma=NULL, varcov=NULL)
{
  localKinhom(X, lambda=lambda, wantL=TRUE, ..., 
              correction=correction, verbose=verbose, rvalue=rvalue,
              sigma=sigma, varcov=varcov)
}

"localK" <-
  function(X, ..., correction="Ripley", verbose=TRUE, rvalue=NULL)
{
  verifyclass(X, "ppp")
  localKengine(X, ..., correction=correction, verbose=verbose, rvalue=rvalue)
}

"localKinhom" <-
  function(X, lambda=NULL, ..., correction="Ripley", verbose=TRUE, rvalue=NULL,
           sigma=NULL, varcov=NULL)
{
  verifyclass(X, "ppp")

  if(is.null(lambda)) {
    # No intensity data provided
    # Estimate density by leave-one-out kernel smoothing
    lambda <- density(X, ..., sigma=sigma, varcov=varcov,
                            at="points", leaveoneout=TRUE)
    lambda <- as.numeric(lambda)
  } else {
    # validate
    if(is.im(lambda)) 
      lambda <- safelookup(lambda, X)
    else if(is.ppm(lambda))
      lambda <- predict(lambda, locations=X, type="trend")
    else if(is.function(lambda)) 
      lambda <- lambda(X$x, X$y)
    else if(is.numeric(lambda) && is.vector(as.numeric(lambda)))
      check.nvector(lambda, npoints(X))
    else stop(paste(sQuote("lambda"),
                    "should be a vector, a pixel image, or a function"))
  }  
  localKengine(X, lambda=lambda, ...,
               correction=correction, verbose=verbose, rvalue=rvalue)
}

"localKengine" <-
  function(X, ..., wantL=FALSE, lambda=NULL,
           correction="Ripley", verbose=TRUE, rvalue=NULL)
{
  npts <- npoints(X)
  W <- X$window
  areaW <- area(W)
  lambda.ave <- npts/areaW
  lambda1.ave <- (npts - 1)/areaW

  weighted <- !is.null(lambda)

  if(is.null(rvalue)) 
    rmaxdefault <- rmax.rule("K", W, lambda.ave)
  else {
    stopifnot(is.numeric(rvalue))
    stopifnot(length(rvalue) == 1)
    stopifnot(rvalue >= 0)
    rmaxdefault <- rvalue
  }
  breaks <- handle.r.b.args(NULL, NULL, W, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max
  
  correction.given <- !missing(correction)
  correction <- pickoption("correction", correction,
                           c(none="none",
                             isotropic="isotropic",
                             Ripley="isotropic",
                             trans="translate",
                             translate="translate",
                             translation="translate",
                             best="best"),
                           multi=FALSE)

  correction <- implemented.for.K(correction, W$type, correction.given)

  # recommended range of r values
  alim <- c(0, min(rmax, rmaxdefault))

  # identify all close pairs
  rmax <- max(r)
  close <- closepairs(X, rmax)
  DIJ <- close$d
  XI <- ppp(close$xi, close$yi, window=W, check=FALSE)
  I <- close$i
  if(weighted) {
    J <- close$j
    lambdaJ <- lambda[J]
    weightJ <- 1/lambdaJ
  } 
  
  # initialise
  df <- as.data.frame(matrix(NA, length(r), npts))
  labl <- desc <- character(npts)

  bkt <- function(x) { paste("[", x, "]", sep="") }

  if(verbose) state <- list()
  
  switch(correction,
         none={
           # uncorrected! For demonstration purposes only!
           for(i in 1:npts) {
             ii <- (I == i)
             wh <- whist(DIJ[ii], breaks$val,
                         if(weighted) weightJ[ii] else NULL)  # no edge weights
             df[,i] <- cumsum(wh)
             icode <- numalign(i, npts)
             names(df)[i] <- paste("un", icode, sep="")
             labl[i] <- paste("%s", bkt(icode), "(r)", sep="")
             desc[i] <- paste("uncorrected estimate of %s",
                              "for point", icode)
             if(verbose) state <- progressreport(i, npts, state=state)
           }
           if(!weighted) df <- df/lambda1.ave
         },
         translate={
           # Translation correction
           XJ <- ppp(close$xj, close$yj, window=W, check=FALSE)
           edgewt <- edge.Trans(XI, XJ, paired=TRUE)
           if(weighted)
             edgewt <- edgewt * weightJ
           for(i in 1:npts) {
             ii <- (I == i)
             wh <- whist(DIJ[ii], breaks$val, edgewt[ii])
             Ktrans <- cumsum(wh)
             df[,i] <- Ktrans
             icode <- numalign(i, npts)
             names(df)[i] <- paste("trans", icode, sep="")
             labl[i] <- paste("%s", bkt(icode), "(r)", sep="")
             desc[i] <- paste("translation-corrected estimate of %s",
                              "for point", icode)
             if(verbose) state <- progressreport(i, npts, state=state)
           }
           if(!weighted) df <- df/lambda1.ave
           h <- diameter(W)/2
           df[r >= h, ] <- NA
         },
         isotropic={
           # Ripley isotropic correction
           edgewt <- edge.Ripley(XI, matrix(DIJ, ncol=1))
           if(weighted)
             edgewt <- edgewt * weightJ
           for(i in 1:npts) {
             ii <- (I == i)
             wh <- whist(DIJ[ii], breaks$val, edgewt[ii])
             Kiso <- cumsum(wh)
             df[,i] <- Kiso
             icode <- numalign(i, npts)
             names(df)[i] <- paste("iso", icode, sep="")
             labl[i] <- paste("%s", bkt(icode), "(r)", sep="")
             desc[i] <- paste("Ripley isotropic correction estimate of %s", 
                              "for point", icode)
             if(verbose) state <- progressreport(i, npts, state=state)
           }
           if(!weighted) df <- df/lambda1.ave
           h <- diameter(W)/2
           df[r >= h, ] <- NA
         })
  # transform values if L required
  if(wantL)
    df <- sqrt(df/pi)
  
  # return vector of values at r=rvalue, if desired
  if(!is.null(rvalue)) {
    nr <- length(r)
    if(r[nr] != rvalue)
      stop("Internal error - rvalue not attained")
    return(as.numeric(df[nr,]))
  }
  # function value table required
  # add r and theo
  if(!wantL) {
    df <- cbind(df, data.frame(r=r, theo=pi * r^2))
    if(!weighted) {
      ylab <- quote(K[loc](r))
      fnam <- "K[loc][',']"
    } else {
      ylab <- quote(Kinhom[loc](r))
      fnam <- "Kinhom[loc][',']"
    }
  } else {
    df <- cbind(df, data.frame(r=r, theo=r))
    if(!weighted) {
      ylab <- quote(L[loc](r))
      fnam <- "L[loc][',']"
    } else {
      ylab <- quote(Linhom[loc](r))
      fnam <- "Linhom[loc][',']"
    }
  }
  desc <- c(desc, c("distance argument r", "theoretical Poisson %s"))
  labl <- c(labl, c("r", "%s[pois](r)"))
  # create fv object
  K <- fv(df, "r", ylab, "theo", , alim, labl, desc, fname=fnam)
  # default is to display them all
  formula(K) <- . ~ r
  unitname(K) <- unitname(X)
  attr(K, "correction") <- correction
  return(K)
}


