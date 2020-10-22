#
#   rlabel.R
#
#   random (re)labelling
#
#   $Revision: 1.12 $   $Date: 2019/04/25 02:59:12 $
#
#
rlabel <- function(X, labels=marks(X), group=NULL, permute=TRUE, nsim=1, drop=TRUE) {
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
    if ( !is.null(group) ) {
      olabels <- order(labels[,group])
      Z <- lapply(Y, FUN= function(x) { 
        ox <- order( marks(x)[,group] )
        sx <- marks(x)[ox,]
        marks(x) <-sx[order(olabels),]
        x
      } )
      Y <- Z
    }
  } else stop("Format of labels argument is not understood")
  return(simulationresult(Y, nsim, drop))
}

