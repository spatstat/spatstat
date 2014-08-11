#
#  residuals.mppm.R
#
# computes residuals for fitted multiple point process model
#
#
#  $Revision: 1.4 $ $Date: 2013/11/11 14:47:19 $
#

residuals.mppm <- function(object, type="raw", ..., 
                          fittedvalues = fitted.mppm(object)) {
  
  verifyclass(object, "mppm")
  userfitted <- !missing(fittedvalues)
  type <- pickoption("type", type,
                     c(inverse="inverse",
                       raw="raw",
                       pearson="pearson",
                       Pearson="pearson"))
  typenames <- c(inverse="inverse-lambda residuals",
                 raw="raw residuals",
                 pearson="Pearson residuals")
  typename <- typenames[[type]]
  
  # Extract quadrature points and weights
  Q <- quad.mppm(object)
  U <- lapply(Q, union.quad) # quadrature point patterns
  Z <- unlist(lapply(Q, is.data)) # indicator data/dummy
  W <- unlist(lapply(Q, w.quad)) # quadrature weights
  # total number of quadrature points
  nquadrature <- length(W)
  # number of quadrature points in each pattern
  nU <- unlist(lapply(Q, n.quad))
  # number of rows of hyperframe
  npat <- object$npat
  # attribution of each quadrature point
  id <- factor(rep(1:npat, nU), levels=1:npat)
  
  # Compute fitted conditional intensity at quadrature points

  if(!is.list(fittedvalues) || length(fittedvalues) != npat)
    stop(paste(sQuote("fittedvalues"), "should be a list of length",
               npat, "containing vectors of fitted values"))
  
  lambda <- unlist(fittedvalues)

  # consistency check
  if(length(lambda) != nquadrature)
    stop(paste(if(!userfitted) "internal error:" else NULL,
               "number of fitted values", paren(length(lambda)),
               "does not match number of quadrature points",
               paren(nquadrature)))

  # indicator is 1 if lambda > 0
  # (adjusted for numerical behaviour of predict.glm)
  indicator <- (lambda > .Machine$double.eps)

  # Evaluate residual measure components
  discrete <- ifelse(Z,
                     switch(type,
                            raw     = 1,
                            inverse = 1/lambda,
                            pearson = 1/sqrt(lambda)
                            ),
                     0)

  density <- switch(type,
                    raw     = -lambda,
                    inverse = -indicator,
                    pearson = -indicator * sqrt(lambda))

  atoms <- as.logical(Z)
  
  # All components
  resdf <- data.frame(discrete=discrete,
                      density=density,
                      atoms=atoms)

  # Split residual data according to point pattern affiliation
  splitres <- split(resdf, id)
  # Associate with quadrature scheme
  reshf <- hyperframe(R=splitres, Q=Q)
  # Convert to signed measures
  answer <- with(reshf, msr(Q, R$discrete[R$atoms], R$density))
  # tag
  answer <- lapply(answer, "attr<-", which="type", value=type)
  answer <- lapply(answer, "attr<-", which="typename", value=typename)
  return(as.listof(answer))
}

