#'
#'    varcount.R
#'
#'   Variance of N(B)
#'
#'  $Revision: 1.8 $  $Date: 2015/11/21 07:02:51 $
#'

varcount <- function(model, B, ..., dimyx=NULL) {
  stopifnot(is.owin(B) || is.im(B) || is.function(B))
  g <- pcfmodel(model)
  if(!is.function(g))
    stop("Pair correlation function cannot be computed")
  if(is.owin(B)) {
    lambdaB <- predict(model, locations=B, ngrid=dimyx, type="intensity")
    v <- varcountEngine(g, B, lambdaB)
  } else {
    f <- if(is.im(B)) B else as.im(B, W=as.owin(model), ..., dimyx=dimyx)
    B <- as.owin(f)
    lambdaB <- predict(model, locations=B, type="intensity")
    v <- varcountEngine(g, B, lambdaB, f)
  } 
  return(v)
}

varcountEngine <- local({

  varcountEngine <- function(g, B, lambdaB, f=1) {
    if(missing(f) || identical(f, 1)) {
      v <- integral(lambdaB) + covterm(g, B, lambdaB)
    } else if(min(f) >= 0) {
      ## nonnegative integrand
      v <- integral(lambdaB * f^2) + covterm(g, B, lambdaB * f)
    } else if(max(f) <= 0) {
      ## nonpositive integrand
      v <- integral(lambdaB * f^2) + covterm(g, B, lambdaB * (-f))
    } else {
      ## integrand has both positive and negative parts
      lamfplus <- eval.im(lambdaB * pmax(0, f))
      lamfminus <- eval.im(lambdaB * pmax(0, -f))
      v <- integral(lambdaB * f^2) +
        (covterm(g, B, lamfplus) + covterm(g, B, lamfminus)
         - covterm(g, B, lamfplus, lamfminus)
         - covterm(g, B, lamfminus, lamfplus))
    }
    return(v)
  }

  covterm <- function(g, B, f, f2) {
    if(missing(f2)) {
      # \int_B \int_B (g(u-v) - 1) f(u) f(v) du dv
      H <- distcdf(B, dW=f)
      a <- integral(f)^2 * (as.numeric(stieltjes(g, H)) - 1)
    } else {
      # \int_B \int_B (g(u-v) - 1) f(u) f2(v) du dv
      H <- distcdf(B, dW=f, dV=f2)
      a <- integral(f) * integral(f2) * (as.numeric(stieltjes(g, H)) - 1)
    }
    return(a)
  }
  
  varcountEngine
})


