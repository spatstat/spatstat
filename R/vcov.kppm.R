#
# vcov.kppm
#
#  vcov method for kppm objects
#
#   Original code: Abdollah Jalilian
#
#   $Revision: 1.1 $  $Date: 2012/02/04 01:42:35 $
#

vcov.kppm <- function(object, ...,
                      what=c("vcov", "corr", "fisher", "internals"))
{
  what <- match.arg(what)
  verifyclass(object, "kppm")
  # extract composite likelihood results
  po <- object$po
  # ensure it was fitted with quadscheme
  if(is.null(getglmfit(po))) {
    warning("Re-fitting model with forcefit=TRUE")
    po <- update(po, forcefit=TRUE)
  }
  # extract quadrature scheme information
  Q <- quad.ppm(po)
  U <- union.quad(Q)
  nU <- npoints(U)
  wt <- w.quad(Q)
  # compute fitted intensity values
  lambda <- fitted(po, type="lambda")
  # extract covariate values
  Z <- model.matrix(po)
  # compute pair correlation function minus 1
  g <- pcfmodel(object)
  r <- pairdist(U)
  gr <- g(r) - 1
  G <- matrix(gr, nU, nU) 
  # evaluate integral
  ff <- Z * lambda * wt
  J <- t(Z) %*% ff
  E <- t(ff) %*% G %*% ff
  # asymptotic covariance matrix in the Poisson case
  J.inv <- try(solve(J))
  # could be singular 
  if(inherits(J.inv, "try-error")) {
    if(what == "internals") {
      return(list(ff=ff, J=J, E=E, J.inv=NULL))
    } else {
      return(NULL)
    }
  }
  # asymptotic covariance matrix in the clustered case
  vc <- J.inv + J.inv %*% E %*% J.inv
  # 
  switch(what,
         vcov={ return(vc) },
         corr={
           sd <- sqrt(diag(vc))
           co <- vc/outer(sd, sd, "*")
           return(co)
         },
         fisher={
           fish <- try(solve(vc))
           if(inherits(fish, "try-error")) fish <- NULL 
           return(fish)
         },
         internals={
           return(list(ff=ff, J=J, E=E, J.inv=J.inv, vc=vc))
         })
  stop(paste("Unrecognised option: what=", what))
}
