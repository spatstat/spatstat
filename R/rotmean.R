##
## rotmean.R
##
## rotational average of pixel values
##
##  $Revision: 1.6 $ $Date: 2014/09/11 07:23:43 $

rotmean <- function(X, ..., origin, Xname, result=c("fv", "im")) {
  if(missing(Xname))
    Xname <- sensiblevarname(short.deparse(substitute(X)), "X")
  trap.extra.arguments(..., .Context="rotmean")
  stopifnot(is.im(X))
  if(!missing(origin))
    X <- shift(X, origin=origin)
  result <- match.arg(result)
  values <- X[drop=TRUE]
  radii <- with(as.data.frame(rasterxy.im(X, drop=TRUE)), sqrt(x^2+y^2))
  ra <- range(radii)
  eps <- sqrt(X$xstep^2 + X$ystep^2)
  a <- unnormdensity(radii,                 from=ra[1], to=ra[2], bw=eps)
  b <- unnormdensity(radii, weights=values, from=ra[1], to=ra[2], bw=eps)
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
  XX <- as.im(X, na.replace=1)
  IM <- as.im(function(x,y,FUN){ FUN(sqrt(x^2+y^2)) }, XX, FUN=FUN)
  return(IM)
}
