##
## rotmean.R
##
## rotational average of pixel values
##
##  $Revision: 1.11 $ $Date: 2018/12/19 05:29:34 $

rotmean <- function(X, ..., origin, padzero=TRUE, Xname, result=c("fv", "im")) {
  if(missing(Xname))
    Xname <- sensiblevarname(short.deparse(substitute(X)), "X")
  trap.extra.arguments(..., .Context="rotmean")
  stopifnot(is.im(X))
  if(!missing(origin)) {
    X <- shift(X, origin=origin)
    backshift <- -getlastshift(X)
  } else {
    backshift <- NULL
  }
  result <- match.arg(result)
  rmax <- with(vertices(Frame(X)), sqrt(max(x^2+y^2)))
  Xunpad <- X
  if(padzero) 
    X <- padimage(na.handle.im(X, 0), 0, W=square(c(-1,1)*rmax))
  Xdata <- as.data.frame(X)
  values <- Xdata$value
  radii <- with(Xdata, sqrt(x^2+y^2))
  ra <- pmin(range(radii), rmax)
  ## eps <- sqrt(X$xstep^2 + X$ystep^2)
  a <- unnormdensity(radii,                 from=ra[1], to=ra[2])
  b <- unnormdensity(radii, weights=values, from=ra[1], to=ra[2], bw=a$bw)
  df <- data.frame(r=a$x, f=b$y/a$y)
  FUN <- fv(df,
            argu="r",
            ylab=substitute(bar(X)(r), list(X=as.name(Xname))),
            valu="f",
            fmla=(. ~ r),
            alim=ra,
            labl=c("r", "%s(r)"),
            desc=c("distance argument r",
                "rotational average"),
            unitname=unitname(X),
            fname=paste0("bar", paren(Xname)))
  attr(FUN, "dotnames") <- "f"
  if(result == "fv") return(FUN)
  ## compute image
  FUN <- as.function(FUN)
  XX <- as.im(Xunpad, na.replace=1)
  IM <- as.im(function(x,y,FUN){ FUN(sqrt(x^2+y^2)) }, XX, FUN=FUN)
  if(!is.null(backshift))
    IM <- shift(IM,backshift)
  return(IM)
}
