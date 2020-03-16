#'  rcelllpp.R
#'
#'  Analogue of Baddeley-Silverman cell process for linear network.
#'
#'  (plus analogue of Switzer's process)
#' 
#'  $Revision: 1.3 $  $Date: 2020/03/16 10:28:51 $

rcelllpp <- local({

  rcelllpp <- function(L, lambda, rnumgen=NULL, ..., saveid=FALSE) {
    if(inherits(L, "lintess")) {
      LT <- L
      L <- as.linnet(LT)
    } else if(inherits(L, "linnet")) {
      #' default tessellation: each segment is a tile
      ns <- nsegments(L)
      df <- data.frame(seg=1:ns, t0=0, t1=1, tile=1:ns)
      LT <- lintess(L, df)
    } else stop("L should be a linnet or lintess")
    #' extract list of tiles
    df <- LT$df
    #' add required data
    df$len <- lengths_psp(as.psp(L))[df$seg]
    #' generate random points
    st <- by(df, df$tile, addpoints, lambda=lambda, rnumgen=rnumgen, ...)
    st <- Reduce(rbind, st)
    X <- lpp(st, L)
    if(saveid)
      attr(X, "cellid") <- marks(cut(X, LT))
    return(X)
  }

  addpoints <- function(df, lambda=1, rnumgen=NULL, ...) {
    #' take a subset of the data frame representing one tile of the tessellation
    #' Add random points in this subset.
    piecelengths <- df$len
    tilelength <- sum(piecelengths)
    mu <- tilelength * lambda
    n <- if(is.null(rnumgen)) rcellnumber(1, mu=mu) else rnumgen(1, mu, ...)
    if(n == 0)
      return(data.frame(seg=integer(0), tp=numeric(0)))
    u <- runif(n, max=tilelength)
    csp <- c(0, cumsum(piecelengths))
    i <- findInterval(u, csp, rightmost.closed=TRUE, all.inside=TRUE)
    seg <- df$seg[i]
    tp <- df$t0[i] + (df$t1 - df$t0)[i] * (u - csp[i])/piecelengths[i]
    return(data.frame(seg=seg, tp=tp))
  }

  rcelllpp
})


rSwitzerlpp <- local({

  rSwitzerlpp <- function(L, lambdacut, rintens=rexp, ...,
                       cuts=c("points", "lines")) {
    stopifnot(inherits(L, "linnet"))
    cuts <- match.arg(cuts)
    switch(cuts,
           points = {
             X <- rpoislpp(lambdacut, L)
             LT <- divide.linnet(X)
           },
           lines = {
             X <- rpoisline(lambdacut, L)
             X <- attr(X, "lines")
             LT <- chop.linnet(L, X)
           })
    Z <- rcelllpp(LT, 1, rNswitzer, rintens=rintens, ...)
    attr(Z, "breaks") <- X
    return(Z)
  }
  
  rNswitzer <- function(n, mu, rintens=rexp, ...) {
    rpois(n, mu * rintens(n, ...))
  }

  rSwitzerlpp
})


  
  
