#'
#' bandwidth selection method of Cronie and Van Lieshout
#' 
#' $Revision: 1.1 $ $Date: 2019/09/30 08:01:20 $

bw.CvL <- function(X, ..., srange=NULL, ns=16, sigma=NULL, warn=TRUE){
  stopifnot(is.ppp(X))
  W <- Window(X)
  areaW <- area.owin(W)
  if(!is.null(sigma)) {
    stopifnot(is.numeric(sigma) && is.vector(sigma))
    ns <- length(sigma)
  } else {
    if(!is.null(srange)) check.range(srange) else {
      nnd <- nndist(X)
      srange <- c(min(nnd[nnd > 0]), diameter(W)/2)
    }
    sigma <- geomseq(from=srange[1L], to=srange[2L], length.out=ns)
  }
  cv <- numeric(ns)
  for(i in 1:ns) {
    si <- sigma[i]
    lamx <- density(X, sigma = si, at = "points",
                    leaveoneout = FALSE, edge = FALSE)
    cv[i] <- ( sum(1/lamx) - areaW )^2
  }
  result <- bw.optim(cv, sigma, iopt=which.min(cv), 
                     creator="bw.CvL",
                     criterion="Cronie and van Lieshout",
                     warnextreme=warn, hargnames="srange",
                     unitname=unitname(X))
  return(result)
}
