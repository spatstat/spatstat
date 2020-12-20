## stienen.R
##
##  Stienen diagram with border correction
##
##  $Revision: 1.9 $ $Date: 2020/12/19 05:25:06 $

stienen <- function(X, ..., bg="grey", border=list(bg=NULL)) {
  Xname <- short.deparse(substitute(X))
  stopifnot(is.ppp(X))
  if(npoints(X) <= 1) {
    W <- Window(X)
    dont.complain.about(W)
    do.call(plot,
            resolve.defaults(list(x=quote(W)),
                             list(...),
                             list(main=Xname)))
    return(invisible(NULL))
  }
  d <- nndist(X)
  b <- bdist.points(X)
  Y <- X %mark% d
  observed <- (b >= d)
  Yobserved <- Y[observed]
  gp <- union(graphicsPars("symbols"), "lwd")
  dont.complain.about(Yobserved)
  do.call.plotfun(plot.ppp,
                  resolve.defaults(list(x=quote(Yobserved),
                                        markscale=1),
                                   list(...),
                                   list(bg=bg),
                                   list(main=Xname)),
                  extrargs=gp)
  if(!identical(border, FALSE)) {
    if(!is.list(border)) border <- list()
    Ycensored <- Y[!observed]
    dont.complain.about(Ycensored)
    do.call.plotfun(plot.ppp,
                    resolve.defaults(list(x=quote(Ycensored),
                                          markscale=1,
                                          add=TRUE),
                                     border,
                                     list(...),
                                     list(bg=bg),
                                     list(cols=grey(0.5), lwd=2)),
                  extrargs=gp)
  }
  return(invisible(NULL))
}

stienenSet <- function(X, edge=TRUE) {
  stopifnot(is.ppp(X))
  nnd <- nndist(X)
  if(!edge) {
    ok <- bdist.points(X) >= nnd
    X <- X[ok]
    nnd <- nnd[ok]
  }
  n <- npoints(X)
  if(n == 0) return(emptywindow(Window(X)))
  if(n == 1) return(Window(X))
  rad <- nnd/2
  if(!all(ok <- (rad > 0))) {
    eps <- min(rad[ok], shortside(Frame(X)))/100
    rad <- pmax(rad, eps)
  }
  delta <- 2 * pi * max(rad)/128
  Z <- disc(rad[1], X[1], delta=delta)
  for(i in 2:n) Z <- union.owin(Z, disc(rad[i], X[i], delta=delta))
  return(Z)
}
