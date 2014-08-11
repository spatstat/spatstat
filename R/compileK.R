# compileK
#
# Function to take a matrix of pairwise distances
# and compile a 'K' function in the format required by spatstat.
#
#   $Revision: 1.5 $  $Date: 2013/04/25 06:37:43 $
# -------------------------------------------------------------------

compileK <- function(D, r, weights=NULL, denom=1, check=TRUE, ratio=FALSE) {
  # process r values
  breaks <- breakpts.from.r(r)
  rmax <- breaks$max
  r    <- breaks$r
  # check that D is a symmetric matrix with nonnegative entries
  if(check)
    stopifnot(is.matrix(D) && isSymmetric(D) && all(D >= 0))
  # ignore the diagonal; throw away any D values greater than rmax
  ok <- (D <= rmax & D > 0)
  Dvalues <- D[ok]
  #
  # weights?
  if(!is.null(weights)) {
    stopifnot(is.matrix(weights) && all(dim(weights)==dim(D)))
    wvalues <- weights[ok]
  } else wvalues <- NULL
  # count the number of D values in each interval (r[k], r[k+1]]
  counts <- whist(Dvalues, breaks=breaks$val, weights=wvalues)
  # cumulative counts: number of D values in [0, r[k])
  Kcount <- cumsum(counts)
  # divide by appropriate denominator
  Kratio <- Kcount/denom
  # wrap it up as an 'fv' object for use in spatstat
  df <- data.frame(r=r, est=Kratio)
  if(!ratio) {
    K <- fv(df, "r", substitute(compileK(r), NULL), "est", . ~ r , c(0,rmax),
            c("r", "%s(r)"),
            c("distance argument r", "estimated %s"),
            fname="compileK")
  } else {
    num <- data.frame(r=r, est=Kcount)
    den <- data.frame(r=r, est=denom)
    K <- ratfv(num, den,
               "r", substitute(compileK(r), NULL), "est", . ~ r , c(0,rmax),
               c("r", "%s(r)"),
               c("distance argument r", "estimated %s"),
               fname="compileK")
  }
  return(K)
}


compilepcf <- function(D, r, weights=NULL, denom=1, check=TRUE,
                       endcorrect=TRUE, ...) {
  # process r values
  breaks <- breakpts.from.r(r)
  if(!breaks$even)
    stop("compilepcf: r values must be evenly spaced", call.=FALSE)
  r    <- breaks$r
  rmax <- breaks$max
  # check that D is a symmetric matrix with nonnegative entries
  if(check)
    stopifnot(is.matrix(D) && isSymmetric(D) && all(D >= 0))
  # ignore the diagonal; throw away any D values greater than rmax
  ok <- (D <= rmax & D > 0)
  Dvalues <- D[ok]
  #
  # weights?
  if(!is.null(weights)) {
    stopifnot(is.matrix(weights) && all(dim(weights)==dim(D)))
    wvalues <- weights[ok]
    totwt <- sum(wvalues)
    normwvalues <- wvalues/totwt
  } else {
    nv <- length(Dvalues)
    normwvalues <- rep.int(1/nv, nv)
    totwt <- nv
  }
  # form kernel estimate
  rmin <- min(r)
  rmax <- max(r)
  nr   <- length(r)
  den <- density(Dvalues, weights=normwvalues,
                 from=rmin, to=rmax, n=nr, ...)
  gval <- den$y * totwt
  # normalise
  gval <- gval/denom
  # edge effect correction at r = 0
  if(endcorrect) {
    one <- do.call("density",
                   resolve.defaults(
                                    list(seq(rmin,rmax,length=512)),
                                    list(bw=den$bw, adjust=1),
                                    list(from=rmin, to=rmax, n=nr),
                                    list(...)))
    onefun <- approxfun(one$x, one$y, rule=2)
    gval <- gval /((rmax-rmin) * onefun(den$x))
  }
  # wrap it up as an 'fv' object for use in spatstat
  df <- data.frame(r=r,
                   est=gval)
  g <- fv(df, "r", substitute(compilepcf(r), NULL), "est", . ~ r , c(0,rmax),
          c("r", "%s(r)"),
          c("distance argument r", "estimated %s"),
          fname="compilepcf")
  attr(g, "bw") <- den$bw
  return(g)
}
