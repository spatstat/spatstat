#'
#'   deldirnet.R
#'
#'   Interfaces to 'deldir' that produce linear networks
#'
#'   $Revision: 1.1 $ $Date: 2020/06/14 10:34:00 $

dirichletNetwork <- function(X, ...) as.linnet(dirichletEdges(X), ...)

delaunayNetwork <- function(X) {
  stopifnot(is.ppp(X))
  X <- unique(X, rule="deldir")
  nX <- npoints(X)
  if(nX == 0) return(NULL)
  if(nX == 1L) return(linnet(X, !diag(TRUE)))
  if(nX == 2L) return(linnet(X, !diag(c(TRUE,TRUE))))
  dd <- safedeldir(X)
  if(is.null(dd)) 
    return(NULL)
  joins <- as.matrix(dd$delsgs[, 5:6])
  return(linnet(X, edges=joins))
}

