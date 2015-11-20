#'
#'    varcount.R
#'
#'   Variance of N(B)
#'
#'  $Revision: 1.6 $  $Date: 2015/11/20 00:47:50 $
#'

varcount <- function(model, B, ..., dimyx=NULL) {
  stopifnot(is.owin(B) || is.im(B) || is.function(B))
  g <- pcfmodel(model)
  if(!is.function(g))
    stop("Pair correlation function cannot be computed")
  if(is.owin(B)) {
    lambdaB <- predict(model, locations=B, ngrid=dimyx, type="intensity")
    varB <- muB <- integral(lambdaB)
    H <- distcdf(B, dW=lambdaB)
    v <- varB + muB^2 * (as.numeric(stieltjes(g, H)) - 1)
  } else {
    f <- if(is.im(B)) B else as.im(B, W=as.owin(model), ..., dimyx=dimyx)
    B <- as.owin(f)
    lambdaB <- predict(model, locations=B, type="intensity")
    if(min(f) >= 0) {
      lamf <- lambdaB * f
      muB <- integral(lamf)
      varB <- integral(lambdaB * f^2)
      H <- distcdf(B, dW=lamf)
      v <- varB + muB^2 * (as.numeric(stieltjes(g, H)) - 1)
    } else if(max(f) <= 0) {
      f <- -f
      lamf <- lambdaB * f
      muB <- integral(lamf)
      varB <- integral(lambdaB * f^2)
      H <- distcdf(B, dW=lamf)
      v <- varB + muB^2 * (as.numeric(stieltjes(g, H)) - 1)
    } else {
      varB <- integral(lambdaB * f^2)
      lamfplus <- lambdaB * eval.im(pmax(0, f))
      lamfminus <- lambdaB * eval.im(pmax(0, -f))
      muBplus <- integral(lamfplus)
      muBminus <- integral(lamfminus)
      Hplusplus <- distcdf(B, dW=lamfplus)
      Hplusminus <- distcdf(B, dW=lamfplus, dV=lamfminus)
      Hminusplus <- distcdf(B, dW=lamfminus, dV=lamfplus)
      Hminusminus <- distcdf(B, dW=lamfminus)
      v <- varB + (
        muBplus^2 * (as.numeric(stieltjes(g, Hplusplus)) - 1)
        + muBminus^2 * (as.numeric(stieltjes(g, Hminusminus)) - 1)
        - muBplus * muBminus * (as.numeric(stieltjes(g, Hplusminus)) - 1)
        - muBminus * muBplus * (as.numeric(stieltjes(g, Hminusplus)) - 1)
        )
    }
  } 
  return(v)
}


