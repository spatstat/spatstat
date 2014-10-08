#
#
#     markcorr.R
#
#     $Revision: 1.65 $ $Date: 2014/10/08 10:47:28 $
#
#    Estimate the mark correlation function
#    and related functions 
#    
# ------------------------------------------------------------------------

"markvario" <-
function(X, correction=c("isotropic", "Ripley", "translate"),
         r=NULL, method="density", ..., normalise=FALSE) {
  m <- onecolumn(marks(X))
  if(!is.numeric(m))
    stop("Marks are not numeric")
  if(missing(correction))
    correction <- NULL
  # Compute estimates
  v <- markcorr(X, f=function(m1, m2) { (1/2) * (m1-m2)^2 },
                r=r, correction=correction, method=method,
                normalise=normalise, ...)
  # adjust theoretical value
  v$theo <- if(normalise) 1 else var(m)    
  # fix labels
  v <- rebadge.fv(v,
                  quote(gamma(r)),
                  "gamma")
  return(v)
}

markconnect <-
function(X, i, j, r=NULL, 
         correction=c("isotropic", "Ripley", "translate"),
         method="density", ..., normalise=FALSE) {
  stopifnot(is.ppp(X) && is.multitype(X))
  if(missing(correction))
    correction <- NULL
  marx <- marks(X)
  lev  <- levels(marx)
  if(missing(i)) i <- lev[1]
  if(missing(j)) j <- lev[2]
  indicateij <- function(m1, m2, i, j) { (m1 == i) & (m2 == j) }
  # compute estimates
  p <- markcorr(X, f=indicateij, r=r,
                correction=correction, method=method,
                ...,
                fargs=list(i=i, j=j),
                normalise=normalise)
  # alter theoretical value and fix labels
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

Emark <- function(X, r=NULL, 
                  correction=c("isotropic", "Ripley", "translate"),
                  method="density", ..., normalise=FALSE) {
  stopifnot(is.ppp(X) && is.marked(X) && is.numeric(marks(X)))
  if(missing(correction))
    correction <- NULL
  f <- function(m1, m2) { m1 }
  E <- markcorr(X, f, r=r,
                correction=correction, method=method,
                ..., normalise=normalise)
  E <- rebadge.fv(E, quote(E(r)), "E")
  return(E)
}

Vmark <- function(X, r=NULL, 
                  correction=c("isotropic", "Ripley", "translate"),
                  method="density", ..., normalise=FALSE) {
  if(missing(correction))
    correction <- NULL
  E <- Emark(X, r=r, correction=correction, method=method, ...,
             normalise=FALSE)
  f2 <- function(m1, m2) { m1^2 }
  E2 <- markcorr(X, f2, r=E$r,
                 correction=correction, method=method,
                 ..., normalise=FALSE)
  V <- eval.fv(E2 - E^2)
  if(normalise) {
    sig2 <- var(marks(X))
    V <- eval.fv(V/sig2)
  }
  V <- rebadge.fv(V, quote(V(r)), "V")
  attr(V, "labl") <- attr(E, "labl")
  return(V)
}

############## workhorses 'markcorr' and 'markcorrint' ####################

markcorrint <-
  function(X, f=NULL, r=NULL, 
           correction=c("isotropic", "Ripley", "translate"), ...,
           f1=NULL, normalise=TRUE, returnL=FALSE, fargs=NULL) {
  # Computes the analogue of Kest(X)
  # where each pair (x_i,x_j) is weighted by w(m_i,m_j)
  #
  # If multiplicative=TRUE then w(u,v) = f(u) f(v)
  # If multiplicative=FALSE then w(u,v) = f(u, v)
  #
  stopifnot(is.ppp(X) && is.marked(X))
  is.marked(X, dfok=FALSE)
  # validate test function
  h <- check.testfun(f, f1, X)
  f     <- h$f
  f1    <- h$f1
  ftype <- h$ftype
  multiplicative <- ftype %in% c("mul", "product")
  # 
  # check corrections
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
  # estimated intensity
  lambda <- X$n/area(Window(X))
  mX <- marks(X)
  switch(ftype,
         mul={
           wt <- mX/lambda
           K <- Kinhom(X, r=r, reciplambda=wt, correction=correction, ...)
           Ef2 <- mean(mX)^2
         },
         equ={
           fXX <- outer(mX, mX, "==")
           wt <- fXX/lambda^2
           K <- Kinhom(X, r=r, reciplambda2=wt, correction=correction, ...)
           mtable <- table(mX)
           Ef2 <- sum(mtable^2)/length(mX)^2
         },
         product={
           f1X <- do.call(f1, append(list(mX), fargs))
           wt <- f1X/lambda
           K <- Kinhom(X, r=r, reciplambda=wt, correction=correction, ...)
           Ef2 <- mean(f1X)^2
         },
         general={
           fXX <- do.call("outer", append(list(mX, mX, f), fargs))
           wt <- fXX/lambda^2
           K <- Kinhom(X, r=r, reciplambda2=wt, correction=correction, ...)
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
           method="density", ..., f1=NULL, normalise=TRUE, fargs=NULL)
{
  # mark correlation function with test function f
  stopifnot(is.ppp(X) && is.marked(X))
  
  # set defaults to NULL
  if(missing(f)) f <- NULL
  if(missing(correction)) correction <- NULL
  
  # handle data frame of marks
  marx <- marks(X, dfok=TRUE)
  if(is.data.frame(marx)) {
    nc <- ncol(marx)
    result <- list()
    for(j in 1:nc) {
      Xj <- X %mark% marx[,j]
      result[[j]] <- markcorr(Xj, f=f, r=r, correction=correction,
                              method=method, ...,
                              f1=f1, normalise=normalise, fargs=fargs)
    }
    result <- as.listof(result)
    names(result) <- colnames(marx)
    return(result)
  }
  
  # validate test function
  h <- check.testfun(f, f1, X)
  f     <- h$f
  f1    <- h$f1
  ftype <- h$ftype
  #
  # 
  npts <- npoints(X)
  W <- X$window

  
  # determine r values 
  rmaxdefault <- rmax.rule("K", W, npts/area(W))
  breaks <- handle.r.b.args(r, NULL, W, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max
        
  if(length(method) > 1)
    stop("Select only one method, please")
  if(method=="density" && !breaks$even)
    stop(paste("Evenly spaced r values are required if method=",
               sQuote("density"), sep=""))
        
  # available selection of edge corrections depends on window
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

  # Denominator
  # Ef = Ef(M,M') when M, M' are independent
  # Apply f to every possible pair of marks, and average
  Ef <- switch(ftype,
               mul = {
                 mean(marx)^2
               },
               equ = {
                 mtable <- table(marx)
                 sum(mtable^2)/sum(mtable)^2
               },
               product={
                 f1m <- do.call(f1, append(list(marx), fargs))
                 mean(f1m)^2
               },
               general = {
                 if(is.null(fargs))
                   mean(outer(marx, marx, f))
                 else
                   mean(do.call("outer", append(list(marx,marx,f),fargs)))
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
    # check validity of denominator
    if(Efdenom == 0)
      stop("Cannot normalise the mark correlation; the denominator is zero")
    else if(Efdenom < 0)
      warning(paste("Problem when normalising the mark correlation:",
                    "the denominator is negative"))
  }
  
  # this will be the output data frame
  result <- data.frame(r=r, theo= rep.int(theory,length(r)))
  desc <- c("distance argument r",
            "theoretical value (independent marks) for %s")
  alim <- c(0, min(rmax, rmaxdefault))
  # determine conventional name of function
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

  # find close pairs of points
  close <- closepairs(X, rmax)
  dIJ <- close$d
  I   <- close$i
  J   <- close$j
  XI <- ppp(close$xi, close$yi, window=W, check=FALSE)

  # apply f to marks of close pairs of points
  #
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

  # check values of f(M1, M2)
  
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
    
  #### Compute estimates ##############
        
  if(any(correction == "translate")) {
    # translation correction
    XJ <- ppp(close$xj, close$yj, window=W, check=FALSE)
    edgewt <- edge.Trans(XI, XJ, paired=TRUE)
    # get smoothed estimate of mark covariance
    Mtrans <- sewsmod(dIJ, ff, edgewt, Efdenom, r, method, ...)
    result <- bind.fv(result,
                      data.frame(trans=Mtrans), "{hat(%s)[%s]^{trans}}(r)",
                      "translation-corrected estimate of %s",
                      "trans")
  }
  if(any(correction == "isotropic")) {
    # Ripley isotropic correction
    edgewt <- edge.Ripley(XI, matrix(dIJ, ncol=1))
    # get smoothed estimate of mark covariance
    Miso <- sewsmod(dIJ, ff, edgewt, Efdenom, r, method, ...)
    result <- bind.fv(result,
                      data.frame(iso=Miso), "{hat(%s)[%s]^{iso}}(r)",
                      "Ripley isotropic correction estimate of %s",
                      "iso")
  }
  # which corrections have been computed?
  nama2 <- names(result)
  corrxns <- rev(nama2[nama2 != "r"])

  # default is to display them all
  formula(result) <- (. ~ r)
  fvnames(result, ".") <- corrxns
  #
  unitname(result) <- unitname(X)
  return(result)
}

sewsmod <- function(d, ff, wt, Ef, rvals, method="smrep", ..., nwtsteps=500) {
  # Smooth Estimate of Weighted Second Moment Density
  # (engine for computing mark correlations, etc)
  # ------
  # Vectors containing one entry for each (close) pair of points
  # d = interpoint distance
  # ff = f(M1, M2) where M1, M2 are marks at the two points
  # wt = edge correction weight
  # -----
  # Ef = E[f(M, M')] where M, M' are independent random marks
  # 
  d <- as.vector(d)
  ff <- as.vector(ff)
  wt <- as.vector(wt)
  switch(method,
         density={
           fw <- ff * wt
           sum.fw <- sum(fw)
           sum.wt <- sum(wt)
           # smooth estimate of kappa_f
           est <- density(d, weights=fw/sum.fw,
                          from=min(rvals), to=max(rvals), n=length(rvals),
                          ...)$y
           numerator <- est * sum.fw
           # smooth estimate of kappa_1
           est0 <- density(d, weights=wt/sum.wt, 
                          from=min(rvals), to=max(rvals), n=length(rvals),
                          ...)$y
           denominator <- est0 * Ef * sum.wt
           result <- numerator/denominator
         },
         sm={
           # This is slow!
           oldopt <- options(warn=-1)
           smok <- require(sm)
           options(oldopt)
           if(!smok)
             stop(paste("Option method=sm requires package sm,",
                        "which is not available"))

           # smooth estimate of kappa_f
           fw <- ff * wt
           est <- sm::sm.density(d, weights=fw,
                                 eval.points=rvals,
                                 display="none", nbins=0, ...)$estimate
           numerator <- est * sum(fw)/sum(est)
           # smooth estimate of kappa_1
           est0 <- sm::sm.density(d, weights=wt,
                                  eval.points=rvals,
                                  display="none", nbins=0, ...)$estimate
           denominator <- est0 * (sum(wt)/ sum(est0)) * Ef
           result <- numerator/denominator
         },
         smrep={
           oldopt <- options(warn=-1)
           smok <- require(sm)
           options(oldopt)
           if(!smok)
             stop(paste("Option method=smrep requires package sm,",
                  "which is not available"))

           hstuff <- resolve.defaults(list(...), list(hmult=1, h.weights=NA))
           if(hstuff$hmult == 1 && all(is.na(hstuff$h.weights)))
             warning("default smoothing parameter may be inappropriate")
           
           # use replication to effect the weights (it's faster)
           nw <- round(nwtsteps * wt/max(wt))
           drep.w <- rep.int(d, nw)
           fw <- ff * wt
           nfw <- round(nwtsteps * fw/max(fw))
           drep.fw <- rep.int(d, nfw)

           # smooth estimate of kappa_f
           est <- sm::sm.density(drep.fw,
                                 eval.points=rvals,
                                 display="none", ...)$estimate
           numerator <- est * sum(fw)/sum(est)
           # smooth estimate of kappa_1
           est0 <- sm::sm.density(drep.w,
                                  eval.points=rvals,
                                  display="none", ...)$estimate
           denominator <- est0 * (sum(wt)/ sum(est0)) * Ef
           result <- numerator/denominator
         },
         loess = {
           # set up data frame
           df <- data.frame(d=d, ff=ff, wt=wt)
           # fit curve to numerator using loess
           fitobj <- loess(ff ~ d, data=df, weights=wt, ...)
           # evaluate fitted curve at desired r values
           Eff <- predict(fitobj, newdata=data.frame(d=rvals))
           # normalise:
           # denominator is the sample mean of all ff[i,j],
           # an estimate of E(ff(M1,M2)) for M1,M2 independent marks
           result <- Eff/Ef
         },
         )
  return(result)
}

############## user interface bits ##################################

check.testfun <- function(f=NULL, f1=NULL, X) {
  # Validate f or f1 as a test function for point pattern X
  # Determine function type 'ftype' ("mul", "equ", "product" or "general")

  fmul <- function(m1, m2) { m1 * m2 }
  fequ <- function(m1, m2) { m1 == m2 }
  f1id <- function(m) { m }

  if(is.null(f) && is.null(f1)) {
    # no functions given
    # default depends on kind of marks
    if(is.multitype(X)) {
      f <- fequ
      ftype <- "equ"
    } else {
      f1 <- f1id
      ftype <- "mul"
    }
  } else if(!is.null(f1)) {
    # f1 given
    # specifies test function of the form f(u,v) = f1(u) f1(v)
    if(!is.null(f))
      warning("argument f ignored (overridden by f1)")
    stopifnot(is.function(f1))
    ftype <- "product"
  } else {
    # f given 
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
      same <- function(f, g) {
        environment(g) <- environment(f)
        identical(f,g)
      }
      ftype <- if(same(f, fmul)) "mul" else if(same(f, fequ)) "equ" else "general"
      if(ftype == "mul" && is.multitype(X))
        stop(paste("Inappropriate choice of function f;",
                   "point pattern is multitype;",
                   "types cannot be multiplied."))
    } else
    stop("Argument f must be a function or the name of a function")
  }
  return(list(f=f, f1=f1, ftype=ftype))
}
