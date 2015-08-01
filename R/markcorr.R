##
##
##     markcorr.R
##
##     $Revision: 1.74 $ $Date: 2015/08/01 08:43:18 $
##
##    Estimate the mark correlation function
##    and related functions 
##    
## ------------------------------------------------------------------------

markvario <- local({

  halfsquarediff <- function(m1, m2) { ((m1-m2)^2)/2 }

  assigntheo <- function(x, value) { x$theo <- value; return(x) }
  
  markvario <- 
    function(X, correction=c("isotropic", "Ripley", "translate"),
             r=NULL, method="density", ..., normalise=FALSE) {
      m <- onecolumn(marks(X))
      if(!is.numeric(m))
        stop("Marks are not numeric")
      if(missing(correction))
        correction <- NULL
      ## Compute estimates
      v <- markcorr(X, f=halfsquarediff, 
                    r=r, correction=correction, method=method,
                    normalise=normalise, ...)
      if(is.fv(v)) v <- anylist(v)
      ## adjust theoretical value and fix labels
      theoval <- if(normalise) 1 else var(m)
      for(i in seq_len(length(v))) {
        v[[i]]$theo <- theoval
        v[[i]] <- rebadge.fv(v[[i]],
                             quote(gamma(r)),
                             "gamma")
      }
      if(length(v) == 1) v <- v[[1]]
      return(v)
    }

  markvario
})

markconnect <- local({

  indicateij <- function(m1, m2, i, j) { (m1 == i) & (m2 == j) }
  
  markconnect <- function(X, i, j, r=NULL, 
                          correction=c("isotropic", "Ripley", "translate"),
                          method="density", ..., normalise=FALSE) {
    stopifnot(is.ppp(X) && is.multitype(X))
    if(missing(correction))
      correction <- NULL
    marx <- marks(X)
    lev  <- levels(marx)
    if(missing(i)) i <- lev[1]
    if(missing(j)) j <- lev[2]
    ## compute estimates
    p <- markcorr(X, f=indicateij, r=r,
                  correction=correction, method=method,
                  ...,
                  fargs=list(i=i, j=j),
                  normalise=normalise)
    ## alter theoretical value and fix labels
    if(!normalise) {
      pipj <- mean(marx==i) * mean(marx==j) 
      p$theo <- pipj
    } else {
      p$theo <- 1
    }
    p <- rebadge.fv(p,
                    new.ylab=substitute(p[i,j](r), list(i=paste(i),j=paste(j))),
                    new.fname=c("p", paste0("list(", i, ",", j, ")")),
                    new.yexp=substitute(p[list(i,j)](r),
                      list(i=paste(i),j=paste(j))))
    return(p)
  }
  markconnect
})


Emark <- local({

  f1 <- function(m1, m2) { m1 }

  Emark <- function(X, r=NULL, 
                    correction=c("isotropic", "Ripley", "translate"),
                    method="density", ..., normalise=FALSE) {
    stopifnot(is.ppp(X) && is.marked(X))
    marx <- marks(X)
    isvec <- is.vector(marx) && is.numeric(marx)
    isdf <- is.data.frame(marx) && all(sapply(as.list(marx), is.numeric))
    if(!(isvec || isdf))
      stop("All marks of X should be numeric")
    if(missing(correction))
      correction <- NULL
    E <- markcorr(X, f1, r=r,
                  correction=correction, method=method,
                  ..., normalise=normalise)
    if(isvec) {
      E <- rebadge.fv(E, quote(E(r)), "E")
    } else {
      E[] <- lapply(E, rebadge.fv, new.ylab=quote(E(r)), new.fname="E")
    }
    return(E)
  }

  Emark
})

Vmark <- local({

  f2 <- function(m1, m2) { m1^2 }

  Vmark <- function(X, r=NULL, 
                    correction=c("isotropic", "Ripley", "translate"),
                    method="density", ..., normalise=FALSE) {
    if(missing(correction))
      correction <- NULL
    E <- Emark(X, r=r, correction=correction, method=method, ...,
             normalise=FALSE)
    E2 <- markcorr(X, f2, r=E$r,
                   correction=correction, method=method,
                   ..., normalise=FALSE)
    if(normalise) 
      sig2 <- var(marks(X))
    if(is.fv(E)) {
      E <- list(E)
      E2 <- list(E2)
    }
    V <- list()
    for(i in seq_along(E)) {
      Ei <- E[[i]]
      E2i <- E2[[i]]
      Vi <- eval.fv(E2i - Ei^2)
      if(normalise) 
        Vi <- eval.fv(Vi/sig2[i,i])
      Vi <- rebadge.fv(Vi, quote(V(r)), "V")
      attr(Vi, "labl") <- attr(Ei, "labl")
      V[[i]] <- Vi
    }
    if(length(V) == 1) return(V[[1]])
    V <- as.anylist(V)
    names(V) <- colnames(marks(X))
    return(V)
  }

  Vmark
})

############## workhorses 'markcorr' and 'markcorrint' ####################

markcorrint <-
Kmark <-
  function(X, f=NULL, r=NULL, 
           correction=c("isotropic", "Ripley", "translate"), ...,
           f1=NULL, normalise=TRUE, returnL=FALSE, fargs=NULL) {
  ## Computes the analogue of Kest(X)
  ## where each pair (x_i,x_j) is weighted by w(m_i,m_j)
  ##
  ## If multiplicative=TRUE then w(u,v) = f(u) f(v)
  ## If multiplicative=FALSE then w(u,v) = f(u, v)
  ##
  stopifnot(is.ppp(X) && is.marked(X))
  is.marked(X, dfok=FALSE)
  W <- Window(X)
  ## 
  if(identical(sys.call()[[1]], as.name('markcorrint')))
    warn.once('markcorrint',
              "markcorrint will be deprecated in future versions of spatstat;",
              "use the equivalent function Kmark")
  ## validate test function
  h <- check.testfun(f, f1, X)
  f     <- h$f
  f1    <- h$f1
  ftype <- h$ftype
  multiplicative <- ftype %in% c("mul", "product")
  ## 
  ## check corrections
  correction.given <- !missing(correction) && !is.null(correction)
  if(is.null(correction))
    correction <- c("isotropic", "Ripley", "translate")
  correction <- pickoption("correction", correction,
                           c(none="none",
                             border="border",
                             "bord.modif"="bord.modif",
                             isotropic="isotropic",
                             Ripley="isotropic",
                             trans="translate",
                             translate="translate",
                             translation="translate",
                             best="best"),
                           multi=TRUE)

  correction <- implemented.for.K(correction, W$type, correction.given)

  isborder  <- correction %in% c("border", "bord.modif")
  if(any(isborder) && !multiplicative) {
    whinge <- paste("Border correction is not valid unless",
                    "test function is of the form f(u,v) = f1(u)*f1(v)")
    correction <- correction[!isborder]
    if(length(correction) == 0)
      stop(whinge)
    else
      warning(whinge)
  }
  ## estimated intensity
  lambda <- intensity(X)
  mX <- marks(X)
  switch(ftype,
         mul={
           wt <- mX/lambda
           K <- Kinhom(X, r=r, reciplambda=wt, correction=correction,
                       ..., renormalise=FALSE)
           Ef2 <- mean(mX)^2
         },
         equ={
           fXX <- outer(mX, mX, "==")
           wt <- fXX/lambda^2
           K <- Kinhom(X, r=r, reciplambda2=wt, correction=correction,
                       ..., renormalise=FALSE)
           mtable <- table(mX)
           Ef2 <- sum(mtable^2)/length(mX)^2
         },
         product={
           f1X <- do.call(f1, append(list(mX), fargs))
           wt <- f1X/lambda
           K <- Kinhom(X, r=r, reciplambda=wt, correction=correction,
                       ..., renormalise=FALSE)
           Ef2 <- mean(f1X)^2
         },
         general={
           fXX <- do.call("outer", append(list(mX, mX, f), fargs))
           wt <- fXX/lambda^2
           K <- Kinhom(X, r=r, reciplambda2=wt, correction=correction,
                       ..., renormalise=FALSE)
           Ef2 <- mean(fXX)
         })
  K$theo <- K$theo * Ef2
  labl <- attr(K, "labl")
  if(normalise)
    K <- eval.fv(K/Ef2)
  if(returnL)
    K <- eval.fv(sqrt(K/pi))
  attr(K, "labl") <- labl
  if(normalise && !returnL) {
    ylab <- quote(K[f](r))
    fnam <- c("K", "f")
  } else if(normalise && returnL) {
    ylab <- quote(L[f](r))
    fnam <- c("L", "f")
  } else if(!normalise && !returnL) {
    ylab <- quote(C[f](r))
    fnam <- c("C", "f")
  } else {
    ylab <- quote(sqrt(C[f](r)/pi))
    fnam <- "sqrt(C[f]/pi)"
  }
  K <- rebadge.fv(K, ylab, fnam)
  return(K)
}

markcorr <-
  function(X, f = function(m1, m2) { m1 * m2}, r=NULL, 
           correction=c("isotropic", "Ripley", "translate"),
           method="density", ...,
           weights=NULL, f1=NULL, normalise=TRUE, fargs=NULL)
{
  ## mark correlation function with test function f
  stopifnot(is.ppp(X) && is.marked(X))
  nX <- npoints(X)
  
  ## set defaults to NULL
  if(missing(f)) f <- NULL
  if(missing(correction)) correction <- NULL

  ## handle data frame of marks
  marx <- marks(X, dfok=TRUE)
  if(is.data.frame(marx)) {
    nc <- ncol(marx)
    result <- list()
    for(j in 1:nc) {
      Xj <- X %mark% marx[,j]
      result[[j]] <- markcorr(Xj, f=f, r=r, correction=correction,
                              method=method, ...,
                              weights=weights,
                              f1=f1, normalise=normalise, fargs=fargs)
    }
    result <- as.anylist(result)
    names(result) <- colnames(marx)
    return(result)
  }
  
  ## weights
  if(unweighted <- is.null(weights)) {
    weights <- rep(1, nX)
  } else {
    stopifnot(is.numeric(weights))
    if(length(weights) == 1) {
      weights <- rep(weights, nX)
    } else check.nvector(weights, nX)
    stopifnot(all(weights > 0))
  }
  
  ## validate test function
  h <- check.testfun(f, f1, X)
  f     <- h$f
  f1    <- h$f1
  ftype <- h$ftype
  ##
  ## 
  npts <- npoints(X)
  W <- X$window
  
  ## determine r values 
  rmaxdefault <- rmax.rule("K", W, npts/area(W))
  breaks <- handle.r.b.args(r, NULL, W, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max
        
  if(length(method) > 1)
    stop("Select only one method, please")
  if(method=="density" && !breaks$even)
    stop(paste("Evenly spaced r values are required if method=",
               sQuote("density"), sep=""))
        
  ## available selection of edge corrections depends on window
  correction.given <- !missing(correction) && !is.null(correction)
  if(is.null(correction))
    correction <- c("isotropic", "Ripley", "translate")
  correction <- pickoption("correction", correction,
                           c(none="none",
                             border="border",
                             "bord.modif"="bord.modif",
                             isotropic="isotropic",
                             Ripley="isotropic",
                             translate="translate",
                             translation="translate",
                             best="best"),
                           multi=TRUE)
  
  correction <- implemented.for.K(correction, W$type, correction.given)

  ## Denominator
  ## Ef = Ef(M,M') when M, M' are independent
  ## Apply f to every possible pair of marks, and average
  Ef <- switch(ftype,
               mul = {
                 mean(marx * weights)^2
               },
               equ = {
                 if(unweighted) {
                   mtable <- table(marx)
                 } else {
                   mtable <- tapply(weights, marx, sum)
                   mtable[is.na(mtable)] <- 0
                 }
                 sum(mtable^2)/nX^2
             },
               product={
                 f1m <- do.call(f1, append(list(marx), fargs))
                 mean(f1m * weights)^2
               },
               general = {
                 mcross <- if(is.null(fargs)) {
                   outer(marx, marx, f)
                 } else {
                   do.call("outer", append(list(marx,marx,f),fargs))
                 }
                 if(unweighted) {
                   mean(mcross)
                 } else {
                   wcross <- outer(weights, weights, "*")
                   mean(mcross * wcross)
                 }
               },
               stop("Internal error: invalid ftype"))

  if(normalise) {
    theory <- 1
    Efdenom <- Ef
  } else {
    theory <- Ef
    Efdenom <- 1
  }

  if(normalise) {
    ## check validity of denominator
    if(Efdenom == 0)
      stop("Cannot normalise the mark correlation; the denominator is zero")
    else if(Efdenom < 0)
      warning(paste("Problem when normalising the mark correlation:",
                    "the denominator is negative"))
  }
  
  ## this will be the output data frame
  result <- data.frame(r=r, theo= rep.int(theory,length(r)))
  desc <- c("distance argument r",
            "theoretical value (independent marks) for %s")
  alim <- c(0, min(rmax, rmaxdefault))
  ## determine conventional name of function
  if(ftype %in% c("mul", "equ")) {
    if(normalise) {
      ylab <- quote(k[mm](r))
      fnam <- c("k", "mm")
    } else {
      ylab <- quote(c[mm](r))
      fnam <- c("c", "mm")
    }
  } else {
    if(normalise) {
      ylab <- quote(k[f](r))
      fnam <- c("k", "f")
    } else {
      ylab <- quote(c[f](r))
      fnam <- c("c", "f")
    }
  }
  result <- fv(result, "r", ylab, "theo", , alim,
               c("r","{%s[%s]^{iid}}(r)"), desc, fname=fnam)

  ## find close pairs of points
  close <- closepairs(X, rmax)
  dIJ <- close$d
  I   <- close$i
  J   <- close$j
  XI <- ppp(close$xi, close$yi, window=W, check=FALSE)

  ## apply f to marks of close pairs of points
  ##
  mI <- marx[I]
  mJ <- marx[J]
  ff <- switch(ftype,
               mul = mI * mJ,
               equ = (mI == mJ),
               product={
                 if(is.null(fargs)) {
                   fI <- f1(mI)
                   fJ <- f1(mJ)
                 } else {
                   fI <- do.call(f1, append(list(mI), fargs))
                   fJ <- do.call(f1, append(list(mJ), fargs))
                 }
                 fI * fJ
               },
               general={
                 if(is.null(fargs))
                   f(marx[I], marx[J])
                 else
                   do.call(f, append(list(marx[I], marx[J]), fargs))
               })

  ## check values of f(M1, M2)
  
  if(is.logical(ff))
    ff <- as.numeric(ff)
  else if(!is.numeric(ff))
    stop("function f did not return numeric values")

  if(any(is.na(ff))) 
    switch(ftype,
           mul=,
           equ=stop("some marks were NA"),
           product=,
           general=stop("function f returned some NA values"))
    
  if(any(ff < 0))
    switch(ftype,
           mul=,
           equ=stop("negative marks are not permitted"),
           product=,
           general=stop("negative values of function f are not permitted"))

  ## weights
  if(!unweighted)
    ff <- ff * weights[I] * weights[J]
  
  #### Compute estimates ##############
        
  if(any(correction == "translate")) {
    ## translation correction
    XJ <- ppp(close$xj, close$yj, window=W, check=FALSE)
    edgewt <- edge.Trans(XI, XJ, paired=TRUE)
    ## get smoothed estimate of mark covariance
    Mtrans <- sewsmod(dIJ, ff, edgewt, Efdenom, r, method, ...)
    result <- bind.fv(result,
                      data.frame(trans=Mtrans), "{hat(%s)[%s]^{trans}}(r)",
                      "translation-corrected estimate of %s",
                      "trans")
  }
  if(any(correction == "isotropic")) {
    ## Ripley isotropic correction
    edgewt <- edge.Ripley(XI, matrix(dIJ, ncol=1))
    ## get smoothed estimate of mark covariance
    Miso <- sewsmod(dIJ, ff, edgewt, Efdenom, r, method, ...)
    result <- bind.fv(result,
                      data.frame(iso=Miso), "{hat(%s)[%s]^{iso}}(r)",
                      "Ripley isotropic correction estimate of %s",
                      "iso")
  }
  ## which corrections have been computed?
  nama2 <- names(result)
  corrxns <- rev(nama2[nama2 != "r"])

  ## default is to display them all
  formula(result) <- (. ~ r)
  fvnames(result, ".") <- corrxns
  ##
  unitname(result) <- unitname(X)
  return(result)
}

## mark cross-correlation function 

markcrosscorr <-
  function(X, r=NULL, 
           correction=c("isotropic", "Ripley", "translate"),
           method="density", ..., normalise=TRUE, Xname=NULL)
{
  if(missing(Xname))
    Xname <- short.deparse(substitute(X))

  stopifnot(is.ppp(X) && is.marked(X))
  npts <- npoints(X)
  W <- Window(X)

  ## available selection of edge corrections depends on window
  correction.given <- !missing(correction) && !is.null(correction)
  if(is.null(correction))
    correction <- c("isotropic", "Ripley", "translate")
  correction <- pickoption("correction", correction,
                           c(none="none",
                             border="border",
                             "bord.modif"="bord.modif",
                             isotropic="isotropic",
                             Ripley="isotropic",
                             translate="translate",
                             translation="translate",
                             best="best"),
                           multi=TRUE)
  correction <- implemented.for.K(correction, W$type, correction.given)
  
  ## determine r values 
  rmaxdefault <- rmax.rule("K", W, npts/area(W))
  breaks <- handle.r.b.args(r, NULL, W, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max
        
  ## find close pairs of points
  close <- closepairs(X, rmax)
  dIJ <- close$d
  I   <- close$i
  J   <- close$j
  XI <- ppp(close$xi, close$yi, window=W, check=FALSE)

  ## determine estimation method
  if(length(method) > 1)
    stop("Select only one method, please")
  if(method=="density" && !breaks$even)
    stop(paste("Evenly spaced r values are required if method=",
               sQuote("density"), sep=""))

  ## ensure marks are a data frame
  marx <- marks(X, dfok=TRUE)
  if(!is.data.frame(marx))
    marx <- data.frame(marks=marx)

  ## convert factor marks to dummy variables
  while(any(isfac <- sapply(marx, is.factor))) {
    i <- min(which(isfac))
    mari <- marx[,i]
    levi <- levels(mari)
    nami <- colnames(marx)[i]
    dumi <- 1 * outer(mari, levi, "==")
    colnames(dumi) <- paste0(nami, levi)
    marx <- as.data.frame(append(marx[,-i,drop=FALSE], list(dumi), after=i-1))
  }
  nc <- ncol(marx)
  nama <- colnames(marx)
  ## loop over all pairs of columns
  funs <- list()
  for(i in 1:nc) {
    marxi <- marx[,i]
    namei <- nama[i]
    for(j in 1:nc) {
      marxj <- marx[,j]
      namej <- nama[j]
      ## Denominator
      ## Ef = E M M' = EM EM'
      ## when M, M' are independent from the respective columns
      Ef <- mean(marxi) * mean(marxj)
      if(normalise) {
        theory <- 1
        Efdenom <- Ef
        ## check validity of denominator
        if(Efdenom == 0)
          stop(paste("Cannot normalise the mark correlation for",
                     namei, "x", namej, "because the denominator is zero"),
               call.=FALSE)
        else if(Efdenom < 0)
          warning(paste("Problem when normalising the mark correlation for",
                        namei, "x", namej,
                        "- the denominator is negative"),
                  call.=FALSE)
      } else {
        theory <- Ef
        Efdenom <- 1
      }
      ## this will be the output data frame
      df.ij <- data.frame(r=r, theo= rep.int(theory,length(r)))
      desc <- c("distance argument r",
                "theoretical value (independent marks) for %s")
      alim <- c(0, min(rmax, rmaxdefault))
      ## determine conventional name of function
      mimj <- as.name(paste0(namei,".",namej))
      if(normalise) {
        ylab <- substitute(k[mm](r), list(mm=mimj))
        fnam <- c("k", as.character(mimj))
      } else {
        ylab <- substitute(c[mm](r), list(mm=mimj))
        fnam <- c("c", as.character(mimj))
      }
      fun.ij <- fv(df.ij, "r", ylab, "theo", , alim,
                   c("r","{%s[%s]^{ind}}(r)"), desc, fname=fnam)

      mI <- marxi[I]
      mJ <- marxj[J]
      ff <- mI * mJ
      ## check values of f(M1, M2)

      if(any(is.na(ff))) 
        stop("some marks were NA", call.=FALSE)

      if(any(ff < 0))
        stop("negative marks are not permitted")
    
      ## Compute estimates ##############
        
      if(any(correction == "translate")) {
        ## translation correction
        XJ <- ppp(close$xj, close$yj, window=W, check=FALSE)
        edgewt <- edge.Trans(XI, XJ, paired=TRUE)
        ## get smoothed estimate of mark covariance
        Mtrans <- sewsmod(dIJ, ff, edgewt, Efdenom, r, method, ...)
        fun.ij <- bind.fv(fun.ij,
                          data.frame(trans=Mtrans),
                          "{hat(%s)[%s]^{trans}}(r)",
                          "translation-corrected estimate of %s",
                          "trans")
      }
      if(any(correction == "isotropic")) {
        ## Ripley isotropic correction
        edgewt <- edge.Ripley(XI, matrix(dIJ, ncol=1))
        ## get smoothed estimate of mark covariance
        Miso <- sewsmod(dIJ, ff, edgewt, Efdenom, r, method, ...)
        fun.ij <- bind.fv(fun.ij,
                          data.frame(iso=Miso), "{hat(%s)[%s]^{iso}}(r)",
                          "Ripley isotropic correction estimate of %s",
                          "iso")
      }
      ## which corrections have been computed?
      nama2 <- names(fun.ij)
      corrxns <- rev(nama2[nama2 != "r"])

      ## default is to display them all
      formula(fun.ij) <- (. ~ r)
      fvnames(fun.ij, ".") <- corrxns
      ##
      unitname(fun.ij) <- unitname(X)
      funs <- append(funs, list(fun.ij))
    }
  }
  # matrix mapping array entries to list positions in 'funs'
  witch <- matrix(1:(nc^2), nc, nc, byrow=TRUE)
  header <- paste("Mark cross-correlation functions for", Xname)
  answer <- fasp(funs, witch, 
                 rowNames=nama, colNames=nama,
                 title=header, dataname=Xname)
  return(answer)
}

sewsmod <- function(d, ff, wt, Ef, rvals, method="smrep", ..., nwtsteps=500) {
  ## Smooth Estimate of Weighted Second Moment Density
  ## (engine for computing mark correlations, etc)
  ## ------
  ## Vectors containing one entry for each (close) pair of points
  ## d = interpoint distance
  ## ff = f(M1, M2) where M1, M2 are marks at the two points
  ## wt = edge correction weight
  ## -----
  ## Ef = E[f(M, M')] where M, M' are independent random marks
  ## 
  d <- as.vector(d)
  ff <- as.vector(ff)
  wt <- as.vector(wt)
  switch(method,
         density={
           fw <- ff * wt
           sum.fw <- sum(fw)
           sum.wt <- sum(wt)
           ## smooth estimate of kappa_f
           est <- density(d, weights=fw/sum.fw,
                          from=min(rvals), to=max(rvals), n=length(rvals),
                          ...)$y
           numerator <- est * sum.fw
           ## smooth estimate of kappa_1
           est0 <- density(d, weights=wt/sum.wt, 
                          from=min(rvals), to=max(rvals), n=length(rvals),
                          ...)$y
           denominator <- est0 * Ef * sum.wt
           result <- numerator/denominator
         },
         sm={
           ## This is slow!
           oldopt <- options(warn=-1)
           smok <- requireNamespace("sm")
           options(oldopt)
           if(!smok)
             stop(paste("Option method=sm requires package sm,",
                        "which is not available"))

           ## smooth estimate of kappa_f
           fw <- ff * wt
           est <- sm::sm.density(d, weights=fw,
                                 eval.points=rvals,
                                 display="none", nbins=0, ...)$estimate
           numerator <- est * sum(fw)/sum(est)
           ## smooth estimate of kappa_1
           est0 <- sm::sm.density(d, weights=wt,
                                  eval.points=rvals,
                                  display="none", nbins=0, ...)$estimate
           denominator <- est0 * (sum(wt)/ sum(est0)) * Ef
           result <- numerator/denominator
         },
         smrep={
           oldopt <- options(warn=-1)
           smok <- requireNamespace("sm")
           options(oldopt)
           if(!smok)
             stop(paste("Option method=smrep requires package sm,",
                  "which is not available"))

           hstuff <- resolve.defaults(list(...), list(hmult=1, h.weights=NA))
           if(hstuff$hmult == 1 && all(is.na(hstuff$h.weights)))
             warning("default smoothing parameter may be inappropriate")
           
           ## use replication to effect the weights (it's faster)
           nw <- round(nwtsteps * wt/max(wt))
           drep.w <- rep.int(d, nw)
           fw <- ff * wt
           nfw <- round(nwtsteps * fw/max(fw))
           drep.fw <- rep.int(d, nfw)

           ## smooth estimate of kappa_f
           est <- sm::sm.density(drep.fw,
                                 eval.points=rvals,
                                 display="none", ...)$estimate
           numerator <- est * sum(fw)/sum(est)
           ## smooth estimate of kappa_1
           est0 <- sm::sm.density(drep.w,
                                  eval.points=rvals,
                                  display="none", ...)$estimate
           denominator <- est0 * (sum(wt)/ sum(est0)) * Ef
           result <- numerator/denominator
         },
         loess = {
           ## set up data frame
           df <- data.frame(d=d, ff=ff, wt=wt)
           ## fit curve to numerator using loess
           fitobj <- loess(ff ~ d, data=df, weights=wt, ...)
           ## evaluate fitted curve at desired r values
           Eff <- predict(fitobj, newdata=data.frame(d=rvals))
           ## normalise:
           ## denominator is the sample mean of all ff[i,j],
           ## an estimate of E(ff(M1,M2)) for M1,M2 independent marks
           result <- Eff/Ef
         },
         )
  return(result)
}

############## user interface bits ##################################

check.testfun <- local({
  
  fmul <- function(m1, m2) { m1 * m2 }
  fequ <- function(m1, m2) { m1 == m2 }
  f1id <- function(m) { m }

  check.testfun <- function(f=NULL, f1=NULL, X) {
    ## Validate f or f1 as a test function for point pattern X
    ## Determine function type 'ftype'
    ##      ("mul", "equ", "product" or "general")

    if(is.null(f) && is.null(f1)) {
      ## no functions given
      ## default depends on kind of marks
      if(is.multitype(X)) {
        f <- fequ
        ftype <- "equ"
      } else {
        f1 <- f1id
        ftype <- "mul"
      }
    } else if(!is.null(f1)) {
      ## f1 given
      ## specifies test function of the form f(u,v) = f1(u) f1(v)
      if(!is.null(f))
        warning("argument f ignored (overridden by f1)")
      stopifnot(is.function(f1))
      ftype <- "product"
    } else {
      ## f given 
      if(is.character(fname <- f)) {
        switch(fname,
               "mul"  = {
                 f1 <- f1id
                 ftype <- "mul"
               },
               "equ" = {
                 f <- fequ
                 ftype <- "equ"
               },
               {
                 f <- get(fname)
                 ftype <- "general"
               })
      } else if(is.function(f)) {
        ftype <- if(isTRUE(all.equal(f, fmul))) "mul" else
                 if(isTRUE(all.equal(f, fequ))) "equ" else "general"
        if(ftype == "mul" && is.multitype(X))
          stop(paste("Inappropriate choice of function f;",
                     "point pattern is multitype;",
                     "types cannot be multiplied."))
      } else
        stop("Argument f must be a function or the name of a function")
    }
    return(list(f=f, f1=f1, ftype=ftype))
  }

  check.testfun
})


