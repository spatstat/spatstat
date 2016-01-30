#'
#'   linnetsurgery.R
#'
#' Surgery on linear networks
#'
#' $Revision: 1.1 $  $Date: 2016/01/30 08:49:00 $
#'

insertVertices <- function(L, ...) {
  L <- as.linnet(L)
  argh <- list(...)
  X <- as.lpp(..., L=L)
  if(npoints(X) == 0) {
    attr(L, "id") <- integer(0)
    return(L)
  }
  co <- coords(X)
  seg <- co$seg
  tp <- co$tp
  v <- L$vertices
  n <- npoints(v)
  nadd <- 0
  vadd <- list(x=numeric(0), y=numeric(0))
  fromadd <- toadd <- id <- integer(0)
  for(theseg in sort(unique(seg))) {
    i <- L$from[theseg]
    j <- L$to[theseg]
    those <- (seg == theseg)
    idthose <- which(those)
    tt <- tp[those]
    oo <- order(tt)
    tt <- tt[oo]
    idadd <- idthose[oo]
    nnew <- length(tt)
    xnew <- with(v, x[i] + tt * diff(x[c(i,j)]))
    ynew <- with(v, y[i] + tt * diff(y[c(i,j)]))
    vnew <- list(x=xnew, y=ynew)
    kk <- n + nadd + (1:nnew)
    fromnew <- c(i, kk)
    tonew   <- c(kk, j)
    nadd <- nadd + nnew
    vadd <- concatxy(vadd, list(x=xnew, y=ynew))
    fromadd <- c(fromadd, fromnew)
    toadd <- c(toadd, tonew)
    id <- c(id, idadd)
  }
  newfrom <- c(L$from[-seg], fromadd)
  newto   <- c(L$to[-seg], toadd)
  newv <- superimpose(v, vadd, check=FALSE)
  Lnew <- linnet(newv, edges=cbind(newfrom, newto),
                 sparse=identical(L$sparse, TRUE))
  newid <- integer(nadd)
  newid[id] <- n + 1:nadd
  attr(Lnew, "id") <- newid
  return(Lnew)
}
