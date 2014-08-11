#
#
#  pcfcross.R
#
#  kernel estimation of cross-type pair correlation function
#
#  Currently computed by differencing pcf
#
pcfcross <- function(X, i, j, ...) {
  stopifnot(is.multitype(X))
  if(missing(i)) i <- levels(marks(X))[1]
  if(missing(j)) j <- levels(marks(X))[2]
  # extract points of types i and j
  Xsplit <- split(X)
  Xi <- Xsplit[[i]]
  Xj <- Xsplit[[j]]
  iname <- make.parseable(paste(i))
  jname <- make.parseable(paste(j))
  if(i == j)  {
    p.ii <- pcf(Xi, ...)
    p.ii <- rebadge.fv(p.ii,
                     new.ylab=substitute(g[i,i](r),
                       list(i=iname)),
                     new.fname=sprintf("g[list(%s,%s)]", iname, iname),
                     new.yexp=substitute(g[list(i,i)](r),
                                      list(i=iname)))
    p.ii <- tweak.fv.entry(p.ii, "theo", new.labl="{%s^{pois}}(r)")
    p.ii <- tweak.fv.entry(p.ii, "trans", new.labl="hat(%s^{trans})(r)")
    p.ii <- tweak.fv.entry(p.ii, "iso", new.labl="hat(%s^{iso})(r)")
    return(p.ii)
  }

  Xall <- superimpose(Xi, Xj, W=X$window)
  # estimate intensities
  lambda.i <- summary(Xi)$intensity
  lambda.j <- summary(Xj)$intensity
  lambda.all <- lambda.i + lambda.j
  # kernel estimates of unmarked pcf's
  p.all <- pcf(Xall, ...)
  rr <- p.all$r
  p.ii   <- do.call("pcf",
                   resolve.defaults(list(Xi),
                                    list(...),
                                    list(r=rr)))
  p.jj   <- do.call("pcf",
                   resolve.defaults(list(Xj),
                                    list(...),
                                    list(r=rr)))
  # differencing
  p.ij <- eval.fv(pmax(0,
                        (p.all * lambda.all^2
                         - p.ii * lambda.i^2
                         - p.jj * lambda.j^2)/(2 * lambda.i * lambda.j)))
  # rebadge
  p.ij <- rebadge.fv(p.ij,
                     new.ylab=substitute(g[i,j](r),
                                         list(i=iname, j=jname)),
                     new.fname=sprintf("g[list(%s,%s)]", iname, jname),
                     new.yexp=substitute(g[list(i,j)](r),
                                         list(i=iname, j=jname)))
  p.ij <- tweak.fv.entry(p.ij, "theo", new.labl="{%s^{pois}}(r)")
  p.ij <- tweak.fv.entry(p.ij, "trans", new.labl="hat(%s^{trans})(r)")
  p.ij <- tweak.fv.entry(p.ij, "iso", new.labl="hat(%s^{iso})(r)")
  return(p.ij)
}

pcfdot <- function(X, i, ...) {
  stopifnot(is.multitype(X))
  marx <- marks(X)
  if(missing(i)) i <- levels(marx)[1]
  # map i and not-i to two types
  marks(X) <- factor(ifelse(marx == i, "i", "n"), levels=c("i", "n"))
  # extract points of type i and not-i
  splitX <- split(X)
  Xi <- splitX[["i"]]
  Xn <- splitX[["n"]]
  Xall <- unmark(X)
  # estimate intensities
  lambda.i <- summary(Xi)$intensity
  lambda.n <- summary(Xn)$intensity
  lambda.all <- lambda.i + lambda.n
  # compute cross type pcf from i to not-i
  p.in <- pcfcross(X, "i", "n", ...)
  rr <- p.in$r
  # compute pcf of type i points using same parameters
  p.ii   <- do.call("pcf",
                   resolve.defaults(list(Xi),
                                    list(...),
                                    list(r=rr)))
  # add
  p.idot <- eval.fv((p.in * lambda.n + p.ii * lambda.i)/lambda.all)
  #
  # rebadge
  iname <- make.parseable(paste(i))
  p.idot <- rebadge.fv(p.idot,
                  substitute(g[i ~ dot](r), list(i=iname)),
                  new.fname=paste("g[", iname, "~ symbol(\"\\267\")]"),
                  new.yexp=substitute(g[i ~ symbol("\267")](r),
                                      list(i=iname)))
  p.idot <- tweak.fv.entry(p.idot, "theo", new.labl="{%s^{pois}}(r)")
  p.idot <- tweak.fv.entry(p.idot, "trans", new.labl="hat(%s^{trans})(r)")
  p.idot <- tweak.fv.entry(p.idot, "iso", new.labl="hat(%s^{iso})(r)")
  return(p.idot)
}

