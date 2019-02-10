#' lurkmppm.R
#'    Lurking variable plot for mppm
#'    $Revision: 1.8 $ $Date: 2019/02/10 08:33:42 $

lurking.mppm <- local({

  zerofun <- function(x) rep(0, length(x))

  threshfun <- function(threshold, value) {
    force(threshold)
    force(value)
    function(x) { value * (x >= threshold) }
  }

  approxcumfun <- function(x, y) {
    stopifnot(length(x) == length(y))
    n <- length(x)
    if(n == 0) return(zerofun)
    if(n == 1) return(threshfun(x, y))
    return(approxfun(x=x, y=y, yleft=0, yright=y[n], rule=2))
  }
  
  as.function.lurk <- function(x, ..., what=c("empirical", "theoretical")) {
    what <- match.arg(what)
    switch(what,
           empirical = {
             with(x$empirical, approxcumfun(covariate, value))
           },
           theoretical = {
             with(x$theoretical, approxcumfun(covariate, mean))
           })
  }

  acceptable <- function(x) { is.im(x) || is.numeric(x) || is.expression(x) }

  approxcumul <- function(yin, xin, xout) {
    if(length(yin) > 1) {
      z <- approx(x=xin, y=yin, xout=xout, rule=2)$y
    } else {
      z <- yin * (xout >= xin)
    }
    return(z)
  }
    
  interpolateworking <- function(object, xx) {
    #' extract working data (variance terms)
    #' and interpolate them at the specified covariate values xx
    w <- attr(object, "working")
    if(is.null(w)) return(NULL)
    w <- as.data.frame(w)
    covariate <- object$theoretical$covariate
    y <- apply(w, 2, approxcumul, xin=covariate, xout=xx)
    return(as.data.frame(y))
  }

  multilurk <- function(object, covariate,
                        type="eem",
                        ...,
                        separate=FALSE,
                        plot.it=TRUE,
                        covname, oldstyle=FALSE, nx=512, main="") {
    cl <- match.call()
    stopifnot(is.mppm(object))
    if(missing(covname)) {
      co <- cl$covariate
      covname <- if(is.name(co)) as.character(co) else
                 if(is.expression(co)) format(co[[1]]) else "covariate"
    }

    Fisher <- vcov(object, what="fisher")
    Vcov <- solve(Fisher)

    if(acceptable(covariate)) {
      cov.is.list <- FALSE
    } else {
      cov.is.list <- is.list(covariate) &&
                     length(covariate) == object$npat &&
                     all(sapply(covariate, acceptable))
      if(!cov.is.list) 
        stop(paste("Argument 'covariate' should be",
                   "a pixel image, a numeric vector, an expression",
                   "or a list of such arguments",
                   "with one entry for each row of original data"),
             call.=FALSE)
    }
    #' pseudo fitted model for each row of data
    futs <- subfits(object)
    #' make lurking variable plot object for each row
    if(cov.is.list) {
      #' list of covariate arguments, one for each row of data
      lurks <- mapply(lurking.ppm,
                      object=futs,
                      covariate=covariate,
                      MoreArgs=list(type=type,
                                    plot.it=FALSE,
                                    ...,
                                    internal=list(saveworking=TRUE,
                                                  Fisher=Fisher),
                                    nx=nx,
                                    oldstyle=oldstyle,
                                    covname=covname),
                      SIMPLIFY=FALSE)
    } else {
      #' One covariate argument to rule them all
      #'    First determine range of covariate values
      covrange <- range(sapply(futs, lurking, covariate=covariate,
                               type=type,
                               internal=list(getrange=TRUE)),
                        na.rm=TRUE)
      #'    Now compute lurking variable plots
      lurks <- anylapply(futs, lurking,
                         covariate=covariate,
                         type=type,
                         plot.it=FALSE, 
                         ...,
                         internal=list(saveworking=TRUE,
                                       Fisher=Fisher,
                                       covrange=covrange),
                         nx=nx,
                         oldstyle=oldstyle,
                         covname=covname)
    }
    if(separate) {
      #' separate lurking variable plots for each row
      if(plot.it) {
        do.call(plot,
                resolve.defaults(list(x=lurks),
                                 list(...),
                                 list(main=main,
                                      mar.panel=c(5,4,2,3))))
        return(invisible(lurks))
      } else {
        return(lurks)
      }
    }
    #' auxiliary info
    infos <- lapply(lurks, attr, which="info")
    #' range of covariate values
    covrange <- range(unlist(lapply(infos, getElement, name="covrange")),
                      na.rm=TRUE)
    xx <- seq(covrange[1], covrange[2], length=nx)
    #' empirical part
    efuns <- lapply(lurks, as.function.lurk, what="empirical")
    vlist <- lapply(efuns, do.call, list(xx))
    sumv <- Reduce("+", vlist)
    empirical <- data.frame(covariate=xx, value=sumv)

    #' similar for theoretical curves
    tfuns <- lapply(lurks, as.function.lurk, what="theoretical")
    vlist <- lapply(tfuns, do.call, list(xx))
    sumv <- Reduce("+", vlist)
    theoretical <- data.frame(covariate=xx, mean=sumv)

    #' variance calculation if available
    wlist <- lapply(lurks, interpolateworking, xx=xx)
    if(!any(sapply(wlist, is.null))) {
      w <- Reduce("+", wlist)
      varI <- w$varI
      if(oldstyle) {
        theoretical$sd <- sqrt(varI)
      } else {
        Bnames <- setdiff(colnames(w), c("varI", "varII"))
        B <- as.matrix(w[, Bnames, drop=FALSE])
        if(ncol(B) != nrow(Vcov)) {
          warning("Internal variance data are incomplete; reverting to oldstyle=TRUE")
          oldstyle <- TRUE
          theoretical$sd <- sqrt(varI)
        } else {
          varII <- quadform(B, Vcov)
          varR <- varI - varII
          ra <- range(varR, finite=TRUE)
          if(ra[1] < 0) {
            warning(paste("Negative values of residual variance!",
                          "Range =",
                          prange(signif(ra, 4))),
                    call.=FALSE)
            varR <- pmax(0, varR)
          }
          theoretical$sd <- sqrt(varR)
        }
      }
    }
    
    ## form result
    result <- list(empirical=empirical, theoretical=theoretical)
    class(result) <- "lurk"

    ## copy info e.g. type of residual
    info <- infos[[1]]
    info$covrange <- covrange
    attr(result, "info") <- info

    if(plot.it)
      plot(result, ..., main=main)
    return(invisible(result))
  }

  multilurk
})

