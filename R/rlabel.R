#
#   rlabel.R
#
#   random (re)labelling
#
#   $Revision: 1.12 $   $Date: 2019/04/25 02:59:12 $
#
#
rlabel <- function(X, labels=marks(X), permute=TRUE, nsim=1, drop=TRUE) {
  stopifnot(is.ppp(X) || is.lpp(X) || is.pp3(X) || is.ppx(X) || is.psp(X))
  if(is.null(labels))
    stop("labels not given and marks not present")
  nthings <- nobjects(X)
  things <- if(is.psp(X)) "segments" else "points"
  if(is.vector(labels) || is.factor(labels)) {
    nlabels <- length(labels)
    if(permute && (nlabels != nthings))
      stop(paste("length of labels vector does not match number of", things))
    Y <- replicate(nsim,
                   X %mark% sample(labels, nthings, replace=!permute),
                   simplify=FALSE)
  } else if(is.data.frame(labels) || is.hyperframe(labels)) {
    nlabels <- nrow(labels)
    if(permute && (nlabels != nthings))
      stop(paste("number of rows of data frame does not match number of",
                 things))
    Y <- replicate(nsim,
                   X %mark% labels[sample(1:nlabels, nthings, replace=!permute), ,drop=FALSE],
                   simplify=FALSE)
  } else stop("Format of labels argument is not understood")
  return(simulationresult(Y, nsim, drop))
}

