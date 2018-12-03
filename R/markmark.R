#'
#'   markmark.R
#'
#'   Mark-mark scatterplot
#'
#'   $Revision: 1.7 $ $Date: 2018/12/03 10:26:38 $


markmarkscatter <- function(X, rmax, ..., col=NULL, symap=NULL, transform=I,
                            jit=FALSE) {
  if(!is.ppp(X) && !is.pp3(X) && !is.ppx(X))
    stop("X should be a point pattern", call.=FALSE)
  if(npoints(X) == 0) {
    warning("Empty point pattern; no plot generated.", call.=FALSE)
    return(invisible(NULL))
  }
  stopifnot(is.marked(X))
  marx <- numeric.columns(marks(X))
  nc <- ncol(marx)
  if(nc == 0)
    stop("No marks are numeric", call.=FALSE)
  if(nc > 1) 
    warning("Multiple columns of numeric marks: using the first column",
            call.=FALSE)
  marx <- marx[,1,drop=TRUE]
  transformed <- !missing(transform)
  marx <- transform(marx)
  if(jit) marx <- jitter(marx, factor=2.5)
  if(is.ppp(X) || is.pp3(X)) {
    cl <- closepairs(X, rmax, what="ijd")
  } else {
    D <- pairdist(X)
    ij <- which(D <= rmax, arr.ind=TRUE)
    cl <- list(i=ij[,1], j=ij[,2], d=as.numeric(D[ij]))
  }
  mi <- marx[cl$i]
  mj <- marx[cl$j]
  d  <- cl$d
  ra <- range(marx)
  Y <- ppp(mi, mj, ra, ra, marks=d, check=FALSE)
  nY <- npoints(Y)
  Y <- Y[order(d, decreasing=TRUE)]
  if(is.null(symap)) {
    if(is.null(col)) 
      col <- grey(seq(0.9, 0, length.out=128))
    if(nY > 0) {
      rd <- c(0, max(d))
      symap <- symbolmap(cols=col, range=rd, size=1, pch=16)
    }
  }
  plot(Y, ..., symap=symap, main="", leg.side="right")
  axis(1)
  axis(2)
  mname <- if(jit && transformed) "Jittered, transformed mark" else
           if(jit) "Jittered mark" else
           if(transformed) "Transformed mark" else "Mark"
  title(xlab=paste(mname, "of first point"),
        ylab=paste(mname, "of second point"))
  if(nY >= 2) {
    mbar2 <- mean(marx)^2
    msd2 <- sqrt(2 * var(marx))
    hyperbola <- function(x) { mbar2/x }
    bandline1 <- function(x) { x + msd2 }
    bandline2 <- function(x) { x - msd2 }
    curve(hyperbola, from=mbar2/ra[2], to=ra[2], add=TRUE)
    curve(bandline1,  from=ra[1], to=ra[2]-msd2, add=TRUE)
    curve(bandline2,  from=ra[1]+msd2, to=ra[2], add=TRUE)
  }
  return(invisible(NULL))
}
