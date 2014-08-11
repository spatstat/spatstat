## percus.R
##
## Percus-Yevick style approximations to pcf and K
##
##  $Revision: 1.4 $ $Date: 2014/01/31 10:10:19 $

pcfmodel.ppm <- local({

  pcfmodel.ppm <- function(model, ...) {
    if(is.multitype(model))
      stop("Not yet implemented for multitype models")
    if(!is.stationary(model))
      stop("Model must be stationary")
    if(is.poisson(model)) return(function(r) rep(1, length(r)))
    inte <- as.interact(model)
    if(inte$family$name != "pairwise")
      stop("Only implemented for pairwise-interaction models")
    lambda <- intensity(model)
    beta <- exp(coef(model)[1])
    par <- inte$par
    pot <- inte$pot
    f <- fitin(model)
    Vcoefs <- f$coefs[f$Vnames]
    Mayer <- inte$Mayer
    G <- Mayer(Vcoefs, inte)
    irange <- reach(inte, epsilon=1e-6)
    G2fun <- inte$Percy
    testit <- resolve.1.default(list(testit=FALSE), list(...))
    if(testit || is.null(G2fun))
      G2fun <- pairwisePercy
    fun <- function(r) {
      pcfapprox(r, beta, lambda, pot, par, Vcoefs, G, G2fun, irange)
    }
    return(fun)
  }

  pcfapprox <- function(r, beta, lambda, pot, par, Vcoefs, G, G2fun, irange) {
    as.numeric((beta/lambda)^2 *
               exp(logpairpot(r, pot, par, Vcoefs)
                   - lambda * G2fun(r, Vcoefs, par, pot=pot,
                                    irange=irange, G=G)))
  }

  logpairpot <- function(r, pot, par, Vcoefs) {
    as.numeric(pot(matrix(r, ncol=1), par) %*% Vcoefs)
  }
  
  negpair <- function(x,y, pot, par, Vcoefs) {
    ## evaluate 1 - g(x,y)
    ## where g(x,y) is pair interaction between (0,0) and (x,y)
    1 - exp(logpairpot(sqrt(x^2+y^2), pot, par, Vcoefs))
  }
  
  pairwisePercy <- function(r, Vcoefs, par, ..., G, pot, irange, dimyx=256) {
    S <- max(max(r), irange)
    ng <- as.im(negpair, square(c(-S,S)),
                  pot=pot, par=par, Vcoefs=Vcoefs,
                  dimyx=dimyx)
    ng2 <- convolve.im(ng)
    rr <- seq(min(r), max(r), length=dimyx[1])
    yy <- ng2[list(x=rr, y=rep.int(0, dimyx[1]))]
    zz <- 2 * G - yy
    z <- approx(rr, zz, r)$y
    return(z)
  }

  pcfmodel.ppm
})

    

Kmodel.ppm <- local({
  
  Kmodel.ppm <- function(model, ...) {
    if(is.poisson(model)) return(function(r) { pi * r^2 })
    pc <- pcfmodel(model, ...)
    K <- function(r) pcf2K(r, pc)
    return(K)
  }

  pcf2K <- function(r, pc) {
    ## integrate the pair correlation function to obtain the K-function
    if(length(r) == 1) {
      ## definite integral
      spcfs <- function(s) { s * pc(s) }
      y <- 2 * pi * integrate(spcfs, lower=0, upper=r)$value
    } else {
      ## indefinite integral
      rr <- seq(0, max(r), length=1025)
      dr <- max(r)/(length(rr) - 1)
      ff <- 2 * pi * rr * pc(rr)
      yy <- dr * cumsum(ff)
      y <- approx(rr, yy, r)$y
    }
    return(y)
  }

  Kmodel.ppm
})
                    
