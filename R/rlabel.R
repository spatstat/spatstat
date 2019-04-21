#
#   rlabel.R
#
#   random (re)labelling
#
#   $Revision: 1.10 $   $Date: 2019/04/21 12:02:38 $
#
#
rlabel <- function(X, labels=marks(X), permute=TRUE, nsim=1, drop=TRUE) {
  stopifnot(is.ppp(X) || is.lpp(X) || is.pp3(X) || is.ppx(X))
  if(is.null(labels))
    stop("labels not given and marks not present")
  npts <- npoints(X)
  if(is.vector(labels) || is.factor(labels)) {
    nlabels <- length(labels)
    if(permute && (nlabels != npts))
      stop("length of labels vector does not match number of points")
    Y <- replicate(nsim,
                   X %mark% sample(labels, npts, replace=!permute),
                   simplify=FALSE)
  } else if(is.data.frame(labels) || is.hyperframe(labels)) {
    nlabels <- nrow(labels)
    if(permute && (nlabels != npts))
      stop("number of rows of data frame does not match number of points")      
    Y <- replicate(nsim,
                   X %mark% labels[sample(1:nlabels, npts, replace=!permute), ,drop=FALSE],
                   simplify=FALSE)
  } else stop("Format of labels argument is not understood")
  if(nsim == 1 && drop)
    return(Y[[1]])
  if(is.sob(X)) Y <- as.solist(Y)
  return(Y)
}

