#'
#'     localKcross.R
#'
#'     original by Ege Rubak
#' 
#'     $Revision: 1.6 $ $Date: 2019/06/23 06:13:31 $

"localLcross" <- function(X, from, to, ..., rmax = NULL, correction = "Ripley") {
  localKcross(X, from, to, ..., rmax = rmax, correction = correction, wantL = TRUE)
}

"localLdot" <- function(X, from, ..., rmax = NULL, correction = "Ripley") {
  localKdot(X, from, ..., rmax = rmax, correction = correction, wantL = TRUE)
}

"localKcross" <-
  function(X, from, to, ..., rmax = NULL, correction="Ripley", verbose=TRUE, rvalue=NULL)
  {
    verifyclass(X, "ppp")
    if(!is.multitype(X, dfok=FALSE)) 
	    stop("Point pattern must be multitype")
    marx <- marks(X)
    if(missing(from))
      from <- levels(marx)[1]
    if(missing(to))
      to <- levels(marx)[2]
    I <- (marx == from)
    if(!any(I))
      stop(paste("No points have mark =", from))
    Iexplain <- paste("points having mark =", from)
    Ikey <- make.parseable(paste(from))
    if(from == to) {
      ## use Kest
      result <- do.call(localK,
                        resolve.defaults(list(X=X[I],
                                              rmax=rmax,
                                              correction=correction,
                                              verbose=verbose,
                                              rvalue=rvalue),
                                         list(...)))
    } else {
      J <- (marx == to)
      if(!any(J))
        stop(paste("No points have mark =", to))
      Jexplain <- paste("points having mark =", to)
      Jkey <- make.parseable(paste(to))
      result <-localKmultiEngine(X, I, J, ...,
                                 Ikey=Ikey, Jkey=Jkey,
                                 Iexplain=Iexplain, Jexplain=Jexplain,
                                 rmax = rmax, correction=correction,
                                 verbose=verbose, rvalue=rvalue)
    }
    return(result)
  }

"localKdot" <- 
function(X, from, ..., rmax = NULL, correction="Ripley", verbose=TRUE, rvalue=NULL)
{
  verifyclass(X, "ppp")
  if(!is.multitype(X, dfok=FALSE)) 
  	stop("Point pattern must be multitype")
  marx <- marks(X)
  if(missing(from)) from <- levels(marx)[1]
  
  I <- (marx == from)
  J <- rep.int(TRUE, X$n)  # i.e. all points
  Iexplain <- paste("points having mark =", from)
  Jexplain <- "points of any type"
  Ikey <- make.parseable(paste(from))
  Jkey <- "."
  
  if(!any(I)) stop(paste("No points have mark =", from))
	
  result <- localKmultiEngine(X, I, J, ...,
                              Iexplain=Iexplain, Jexplain=Jexplain,
                              Ikey=Ikey, Jkey=Jkey,
                              rmax = rmax, correction=correction,
                              verbose=verbose, rvalue=rvalue)
  attr(result, "indices") <- list(from=from)
  return(result)
}

"localKcross.inhom" <-
  function(X, from, to, lambdaFrom=NULL, lambdaTo=NULL, ..., rmax = NULL,
           correction = "Ripley",
           sigma=NULL, varcov=NULL,
           lambdaX=NULL, update=TRUE, leaveoneout=TRUE)
  {
    verifyclass(X, "ppp")
    if(!is.multitype(X, dfok=FALSE))
      stop("Point pattern must be multitype")
    miss.update <- missing(update)
    miss.leave <- missing(leaveoneout)
    marx <- marks(X)
    if(missing(from))
      from <- levels(marx)[1]
    if(missing(to))
      to <- levels(marx)[2]
    I <- (marx == from)
    J <- (marx == to)
    Iexplain <- paste("points having mark =", from)
    Jexplain <- paste("points having mark =", to)
    Ikey <- make.parseable(paste(from))
    Jkey <- make.parseable(paste(to))
    K <- localKmultiEngine(X, I, J, lambdaFrom, lambdaTo, ..., rmax = rmax,
                           Iexplain=Iexplain, Jexplain=Jexplain,
                           Ikey=Ikey, Jkey=Jkey,
                           correction=correction,
                           sigma=sigma, varcov=varcov,
                           lambdaX=lambdaX, update=update,
                           leaveoneout=leaveoneout,
                           miss.update=miss.update, miss.leave=miss.leave)
    attr(K, "indices") <- list(from=from, to=to)
    return(K)
  }

localLcross.inhom <- function(X, from, to, lambdaFrom = NULL, lambdaTo = NULL, ..., rmax = NULL) {
  localKcross.inhom(X, from, to, lambdaFrom, lambdaTo, ..., rmax = rmax, wantL = TRUE)
}

"localKmultiEngine" <-
  function(X, from, to, lambdaFrom=NULL, lambdaTo=NULL, ...,
           rmax = NULL, wantL=FALSE,
           correction="Ripley", verbose=TRUE, rvalue=NULL,
           sigma=NULL, varcov=NULL,
           lambdaX=NULL, update=TRUE, leaveoneout=TRUE,
           Iexplain="points satisfying condition I",
           Jexplain="points satisfying condition J",
           Ikey="I",
           Jkey="J",
           miss.update=missing(update),
           miss.leave=missing(leaveoneout))
  {
    npts <- npoints(X)
    W <- Window(X)
    areaW <- area(W)
    lambda.ave <- npts/areaW
    
    from <- ppsubset(X, from)
    to <- ppsubset(X, to)
    if(is.null(from) || is.null(to))
      stop("from and to must be valid subset indices")
    
    if(!any(from)) stop("no points belong to subset from")
    if(!any(to)) stop("no points belong to subset to")
    
    X_from <- X[from]
    X_to <- X[to]
    
    n_from <- sum(from)
    n_to <- sum(to)
    
    lambdaFrom.ave <- n_from/areaW
    lambdaTo.ave <- n_to/areaW
    
    weighted <- !is.null(lambdaFrom)
    if(weighted){
      stopifnot(!is.null(lambdaTo))
      lambdas <- resolve.lambda.cross(X, from, to, lambdaFrom, lambdaTo, ...,
                                      lambdaX = lambdaX,
                                      sigma = sigma, varcov = varcov,
                                      leaveoneout = leaveoneout,
                                      update = update,
                                      Iexplain=Iexplain,
                                      Jexplain=Jexplain,
                                      miss.update=miss.update,
                                      miss.leave=miss.leave,
                                      caller = "localKcrossEngine")
      lambdaFrom <- lambdas$lambdaI
      lambdaTo <- lambdas$lambdaJ
    }
    
    if(is.null(rvalue)) 
      rmaxdefault <- rmax %orifnull% rmax.rule("K", W, lambda.ave)
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
    close <- crosspairs(X_from, X_to, rmax)
    # close$i and close$j are serial numbers in X_from and X_to respectively;        
    # map them to original serial numbers in X
    orig <- seq_len(npts)
    imap <- orig[from]
    jmap <- orig[to]
    iX <- imap[close$i]
    jX <- jmap[close$j]
    # eliminate any identical pairs
    if(any(from & to)) {
      ok <- (iX != jX)
      if(!all(ok)) {
        close$i  <- close$i[ok]
        close$j  <- close$j[ok]
        close$d  <- close$d[ok]
        close$xi  <- close$xi[ok]
        close$xj  <- close$xj[ok]
        close$yi  <- close$yi[ok]
        close$yj  <- close$yj[ok]
      }
    }
    # extract information for these pairs (relative to orderings of X_from, X_to)
    DIJ <- close$d
    XI <- ppp(close$xi, close$yi, window=W, check=FALSE)
    I <- close$i
    J <- close$j
    if(weighted) {
      lambdaI <- lambdaFrom[I]
      lambdaJ <- lambdaTo[J]
      weightI <- 1/lambdaI
      weightJ <- 1/lambdaJ
    } 
    
    # initialise
    df <- as.data.frame(matrix(NA, length(r), n_from))
    labl <- desc <- character(n_from)
    
    if(verbose) state <- list()
    
    switch(correction,
           none={
             # uncorrected! For demonstration purposes only!
             for(i in 1:n_from) {
               ii <- (I == i)
               ## Below
               wh <- whist(DIJ[ii], breaks$val,
                           if(weighted) weightJ[ii] else NULL)  # no edge weights
               Knone <- cumsum(wh)
               ## Tweaking factor to express Kcross.inhom as unweighted average of local contrib.
               if(weighted) Knone <- Knone * lambdaFrom.ave/lambdaFrom[i]
               df[,i] <- Knone
               icode <- numalign(i, n_from)
               names(df)[i] <- paste("un", icode, sep="")
               labl[i] <- makefvlabel(NULL, "hat", character(2), icode)
               desc[i] <- paste("uncorrected estimate of %s",
                                "for point", icode)
               if(verbose) state <- progressreport(i, n_from, state=state)
               
             }
             if(!weighted) df <- df/lambdaTo.ave
           },
           translate={
             # Translation correction
             XJ <- ppp(close$xj, close$yj, window=W, check=FALSE)
             edgewt <- edge.Trans(XI, XJ, paired=TRUE)
             if(weighted)
               edgewt <- edgewt * weightJ
             for(i in 1:n_from) {
               ii <- (I == i)
               wh <- whist(DIJ[ii], breaks$val, edgewt[ii])
               Ktrans <- cumsum(wh)
               ## Tweaking factor to express Kcross.inhom as unweighted average of local contrib.
               if(weighted) Ktrans <- Ktrans * lambdaFrom.ave/lambdaFrom[i]
               df[,i] <- Ktrans
               icode <- numalign(i, n_from)
               names(df)[i] <- paste("trans", icode, sep="")
               labl[i] <- makefvlabel(NULL, "hat", character(2), icode)
               desc[i] <- paste("translation-corrected estimate of %s",
                                "for point", icode)
               if(verbose) state <- progressreport(i, n_from, state=state)
             }
             if(!weighted) df <- df/lambdaTo.ave
             h <- diameter(W)/2
             df[r >= h, ] <- NA
           },
           isotropic={
             # Ripley isotropic correction
             edgewt <- edge.Ripley(XI, matrix(DIJ, ncol=1))
             if(weighted)
               edgewt <- edgewt * weightJ
             for(i in 1:n_from) {
               ii <- (I == i)
               wh <- whist(DIJ[ii], breaks$val, edgewt[ii])
               Kiso <- cumsum(wh)
               ## Tweaking factor to express Kcross.inhom as unweighted average of local contrib.
               if(weighted) Kiso <- Kiso * lambdaFrom.ave/lambdaFrom[i]
               df[,i] <- Kiso
               icode <- numalign(i, n_from)
               names(df)[i] <- paste("iso", icode, sep="")
               labl[i] <- makefvlabel(NULL, "hat", character(2), icode)
               desc[i] <- paste("Ripley isotropic correction estimate of %s", 
                                "for point", icode)
               if(verbose) state <- progressreport(i, n_from, state=state)
             }
             if(!weighted) df <- df/lambdaTo.ave
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
    ## function value table required
    ## add r and theo
    df <- cbind(df,
                data.frame(r=r,
                           theo=if(wantL) r else (pi * r^2)))
    desc <- c(desc, c("distance argument r", "theoretical Poisson %s"))
    labl <- c(labl, c("r", "{%s[%s]^{pois}}(r)"))
    ## Handle 'dot' symbol
    if(identical(Jkey, ".")) {
      Jkeyname <- "symbol(\"\\267\")"
      Jkeylab  <- quote(dot)
      Jkeyexpr <- quote(symbol("\267"))
    } else Jkeyname <- Jkeylab <- Jkeyexpr <- Jkey
    ## Determine fv labels
    if(!wantL) {
      if(!weighted) {
        fnam <- c("K", paste0("list(loc,", Ikey, ",", Jkeyname, ")"))
        ylab <- substitute(K[loc,I,J](r), list(I=Ikey, J=Jkeylab))
        yexp <- substitute(K[list(loc,I,J)](r), list(I=Ikey, J=Jkeyexpr))
      } else {
        fnam <- c("K", paste0("list(inhom,loc,", Ikey, ",", Jkeyname, ")"))
        ylab <- substitute(K[inhom,loc,I,J](r), list(I=Ikey, J=Jkeylab))
        yexp <- substitute(K[list(inhom,loc,I,J)](r), list(I=Ikey, J=Jkeyexpr))
      }
    } else {
      if(!weighted) {
        fnam <- c("L", paste0("list(loc,", Ikey, ",", Jkeyname, ")"))
        ylab <- substitute(L[loc,I,J](r), list(I=Ikey, J=Jkeylab))
        yexp <- substitute(L[list(loc,I,J)](r), list(I=Ikey, J=Jkeyexpr))
      } else {
        fnam <- c("L", paste0("list(inhom,loc,", Ikey, ",", Jkeyname, ")"))
        ylab <- substitute(L[inhom,loc,I,J](r), list(I=Ikey, J=Jkeylab))
        yexp <- substitute(L[list(inhom,loc,I,J)](r), list(I=Ikey, J=Jkeyexpr))
      }
    }
    # create fv object
    K <- fv(df, "r", ylab, "theo", , alim, labl, desc,
            fname=fnam, yexp=yexp)
    # default is to display them all
    formula(K) <- . ~ r
    unitname(K) <- unitname(X)
    attr(K, "correction") <- correction
    if(weighted && lambdas$danger)
      attr(K, "dangerous") <- lambdas$dangerous
    ### TEMPORARY HACK to save info about the "from" points
    attr(K, "Xfrom") <- X_from
    return(K)
  }

resolve.lambda.cross <- function(X, I, J, lambdaI, lambdaJ, ..., lambdaX,
                                 sigma, varcov, leaveoneout, update,
                                 lambdaIJ=NULL,
                                 Iexplain="points satisfying condition I",
                                 Jexplain="points satisfying condition J",
                                 miss.update=missing(update),
                                 miss.leave=missing(leaveoneout),
                                 caller){
  dangerous <- c("lambdaI", "lambdaJ")
  dangerI <- dangerJ <- TRUE
  XI <- X[I]
  XJ <- X[J]
  nI <- npoints(XI)
  nJ <- npoints(XJ)

  if(!is.null(lambdaX)) {
    ## Intensity values for all points of X
    if(!is.null(lambdaI))
      warning("lambdaI was ignored, because lambdaX was given", call.=FALSE)
    if(!is.null(lambdaJ))
      warning("lambdaJ was ignored, because lambdaX was given", call.=FALSE)
    if(is.im(lambdaX)) {
      ## Look up intensity values
      lambdaI <- safelookup(lambdaX, XI)
      lambdaJ <- safelookup(lambdaX, XJ)
    } else if(is.function(lambdaX)) {
      ## evaluate function at locations
      lambdaI <- lambdaX(XI$x, XI$y)
      lambdaJ <- lambdaX(XJ$x, XJ$y)
    } else if(is.numeric(lambdaX) && is.vector(as.numeric(lambdaX))) {
      ## vector of intensity values
      if(length(lambdaX) != npoints(X))
        stop(paste("The length of", sQuote("lambdaX"),
                   "should equal the number of points of X"))
      lambdaI <- lambdaX[I]
      lambdaJ <- lambdaX[J]
    } else if(is.ppm(lambdaX) || is.kppm(lambdaX) || is.dppm(lambdaX)) {
      ## point process model provides intensity
      model <- lambdaX
      if(!update) {
        ## just use intensity of fitted model
        lambdaI <- predict(model, locations=XI, type="trend")
        lambdaJ <- predict(model, locations=XJ, type="trend")
      } else {
        ## re-fit model to data X
        if(is.ppm(model)) {
          model <- update(model, Q=X)
          lambdaX <- fitted(model, dataonly=TRUE, leaveoneout=leaveoneout)
        } else {
          model <- update(model, X=X)
          lambdaX <- fitted(model, dataonly=TRUE, leaveoneout=leaveoneout)
        }
        lambdaI <- lambdaX[I]
        lambdaJ <- lambdaX[J]
        dangerI <- dangerJ <- FALSE
        dangerous <- "lambdaIJ"
        if(miss.update & caller == "Kmulti.inhom") 
          warn.once(key="Kmulti.inhom.update",
                    "The behaviour of Kmulti.inhom when lambda is a ppm object",
                    "has changed (in spatstat 1.45-3 and later).",
                    "See help(Kmulti.inhom)")
      }
    } else stop(paste("Argument lambdaX is not understood:",
                      "it should be a numeric vector,",
                      "an image, a function(x,y)",
                      "or a fitted point process model (ppm, kppm or dppm)"))
  } else {
    ## lambdaI, lambdaJ expected
    if(is.null(lambdaI)) {
      ## estimate intensity
      dangerI <- FALSE
      dangerous <- setdiff(dangerous, "lambdaI")
      lambdaI <- density(XI, ..., sigma=sigma, varcov=varcov,
                         at="points", leaveoneout=leaveoneout)
    } else if(is.im(lambdaI)) {
      ## look up intensity values
      lambdaI <- safelookup(lambdaI, XI)
    } else if(is.function(lambdaI)) {
      ## evaluate function at locations
      lambdaI <- lambdaI(XI$x, XI$y)
    } else if(is.numeric(lambdaI) && is.vector(as.numeric(lambdaI))) {
      ## validate intensity vector
      check.nvector(lambdaI, nI, things=Iexplain)
    } else if(is.ppm(lambdaI) || is.kppm(lambdaI) || is.dppm(lambdaI)) {
      ## point process model provides intensity
      model <- lambdaI
      if(!update) {
        ## just use intensity of fitted model
        lambdaI <- predict(model, locations=XI, type="trend")
      } else {
        ## re-fit model to data X
        model <- if(is.ppm(model)) update(model, Q=X) else update(model, X=X)
        lambdaX <- fitted(model, dataonly=TRUE, leaveoneout=leaveoneout)
        lambdaI <- lambdaX[I]
        dangerI <- FALSE
        dangerous <- setdiff(dangerous, "lambdaI")
        if(miss.update && caller == "Kmulti.inhom")
          warn.once(key="Kmulti.inhom.update",
                    "The behaviour of Kmulti.inhom when lambda is a ppm object",
                    "has changed (in spatstat 1.45-3 and later).",
                    "See help(Kmulti.inhom)")
      }
    } else stop(paste(sQuote("lambdaI"), "should be a vector or an image"))
    
    if(is.null(lambdaJ)) {
      ## estimate intensity
      dangerJ <- FALSE
      dangerous <- setdiff(dangerous, "lambdaJ")
      lambdaJ <- density(XJ, ..., sigma=sigma, varcov=varcov,
                         at="points", leaveoneout=leaveoneout)
    } else if(is.im(lambdaJ)) {
      ## look up intensity values
      lambdaJ <- safelookup(lambdaJ, XJ)
    } else if(is.function(lambdaJ)) {
      ## evaluate function at locations
      lambdaJ <- lambdaJ(XJ$x, XJ$y)
    } else if(is.numeric(lambdaJ) && is.vector(as.numeric(lambdaJ))) {
      ## validate intensity vector
      check.nvector(lambdaJ, nJ, things=Jexplain)
    } else if(is.ppm(lambdaJ) || is.kppm(lambdaJ) || is.dppm(lambdaJ)) {
      ## point process model provides intensity
      model <- lambdaJ
      if(!update) {
        ## just use intensity of fitted model
        lambdaJ <- predict(model, locations=XJ, type="trend")
      } else {
        ## re-fit model to data X
        model <- if(is.ppm(model)) update(model, Q=X) else update(model, X=X)
        lambdaX <- fitted(model, dataonly=TRUE, leaveoneout=leaveoneout)
        lambdaJ <- lambdaX[J]
        dangerJ <- FALSE
        dangerous <- setdiff(dangerous, "lambdaJ")
        if(miss.update & caller == "Kmulti.inhom")
          warn.once(key="Kmulti.inhom.update",
                    "The behaviour of Kmulti.inhom when lambda is a ppm object",
                    "has changed (in spatstat 1.45-3 and later).",
                    "See help(Kmulti.inhom)")
      }
    } else 
      stop(paste(sQuote("lambdaJ"), "should be a vector or an image"))
    
    ## Weight for each pair
    if(!is.null(lambdaIJ)) {
      dangerIJ <- TRUE
      dangerous <- union(dangerous, "lambdaIJ")
      if(!is.matrix(lambdaIJ))
        stop("lambdaIJ should be a matrix")
      if(nrow(lambdaIJ) != nI)
        stop(paste("nrow(lambdaIJ) should equal the number of", Iexplain))
      if(ncol(lambdaIJ) != nJ)
        stop(paste("ncol(lambdaIJ) should equal the number of", Jexplain))
    } else {
      dangerIJ <- FALSE
    }
    
    danger <- dangerI || dangerJ || dangerIJ
    
    return(list(lambdaI = lambdaI, lambdaJ = lambdaJ, lambdaIJ=lambdaIJ,
                danger = danger, dangerous = dangerous))
  }
}
