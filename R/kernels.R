#
# kernels.R
#
#  rXXX, dXXX, pXXX and qXXX for kernels
#
#  $Revision: 1.19 $  $Date: 2018/06/07 05:42:54 $
#

match.kernel <- function(kernel) {
  kernel.map <- c(Gaussian    ="gaussian",
                  gaussian    ="gaussian",
                  Normal      ="gaussian",
                  normal      ="gaussian",
                  rectangular ="rectangular",
                  triangular  ="triangular",
                  Epanechnikov="epanechnikov",
                  epanechnikov="epanechnikov",
                  biweight    ="biweight",
                  cosine      ="cosine",
                  optcosine   ="optcosine"
                  )
  ker <- pickoption("kernel", kernel, kernel.map)
  return(ker)
}

kernel.factor <- function(kernel="gaussian") {
  # This function returns the factor c such that
  #              h = c * sigma
  # where sigma is the standard deviation of the kernel, and
  # h is the corresponding bandwidth parameter as conventionally defined.

  # Conventionally h is defined as a scale factor
  # relative to the `standard form' of the kernel, namely the 
  # form with support [-1,1], except in the Gaussian case where
  # the standard form is N(0,1).
  
  # Thus the standard form of the kernel (h=1) has standard deviation 1/c.
  
  # The kernel with standard deviation 1 has support [-c,c]
  # except for gaussian case.
  
  kernel <- match.kernel(kernel)
  switch(kernel,
         gaussian     = 1,
         rectangular  = sqrt(3),
         triangular   = sqrt(6),
         epanechnikov = sqrt(5),
         biweight     = sqrt(7),
         cosine       = 1/sqrt(1/3 - 2/pi^2),
         optcosine    = 1/sqrt(1 - 8/pi^2))
}

rkernel <- function(n, kernel="gaussian", mean=0, sd=1) {
  kernel <- match.kernel(kernel)
  if(kernel == "gaussian")
    return(rnorm(n, mean=mean, sd=sd))
  # inverse cdf transformation
  u <- runif(n)
  qkernel(u, kernel, mean=mean, sd=sd)
}

dkernel <- function(x,  kernel="gaussian", mean=0, sd=1) {
  kernel <- match.kernel(kernel)
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(sd) && length(sd) == 1 && sd > 0)
  a <- sd * kernel.factor(kernel)
  y <- abs(x-mean)/a
  dens <-
    switch(kernel,
           gaussian       = { dnorm(y) },
           rectangular    = { ifelse(y < 1, 1/2, 0) },
           triangular     = { ifelse(y < 1, (1 - y), 0) },
           epanechnikov   = { ifelse(y < 1, (3/4) * (1 - y^2), 0) },
           biweight       = { ifelse(y < 1, (15/16) * (1 - y^2)^2, 0) },
           cosine         = { ifelse(y < 1, (1 + cos(pi * y))/2, 0) },
           optcosine      = { ifelse(y < 1, (pi/4) * cos(pi * y/2), 0) }
           )
  dens/a
}

pkernel <- function(q, kernel="gaussian", mean=0, sd=1, lower.tail=TRUE){
  kernel <- match.kernel(kernel)
  stopifnot(is.numeric(q))
  stopifnot(is.numeric(sd) && length(sd) == 1 && sd > 0)
  a <- sd * kernel.factor(kernel)
  y <- (q-mean)/a
  switch(kernel,
         gaussian = {
           pnorm(y, lower.tail=lower.tail)
         },
         rectangular = {
           punif(y, min=-1, max=1, lower.tail=lower.tail)
         },
         triangular = {
           p <- ifelse(y < -1, 0, ifelse(y > 1, 1,
                       ifelse(y < 0, y + y^2/2 + 1/2,
                              y - y^2/2 + 1/2)))
           if(lower.tail) p else (1 - p)
         },
         epanechnikov = {
           p <- ifelse(y < -1, 0, ifelse(y > 1, 1,
                        (2 + 3 * y - y^3)/4))
           if(lower.tail) p else (1 - p)
         },
         biweight = {
           p <- ifelse(y < -1, 0, ifelse(y > 1, 1,
                       (15 * y - 10 * y^3 + 3 * y^5 + 8)/16))
           if(lower.tail) p else (1 - p)
         },
         cosine = {
           p <- ifelse(y < -1, 0, ifelse(y > 1, 1,
                       (y + sin(pi * y)/pi + 1)/2))
           if(lower.tail) p else (1 - p)
         },
         optcosine = {
           p <- ifelse(y < -1, 0, ifelse(y > 1, 1,
                       (sin(pi * y/2) + 1)/2))
           if(lower.tail) p else (1 - p)
         })
}

qkernel <- function(p, kernel="gaussian", mean=0, sd=1, lower.tail=TRUE) {
  kernel <- match.kernel(kernel)
  stopifnot(is.numeric(p))
  stopifnot(is.numeric(sd) && length(sd) == 1 && sd > 0)
  a <- sd * kernel.factor(kernel)
  if(!lower.tail)
    p <- 1 - p
  y <-
    switch(kernel,
           gaussian = {
             qnorm(p, lower.tail=lower.tail)
           },
           rectangular = {
             qunif(p, min=-1, max=1, lower.tail=lower.tail)
           },
           triangular = {
             ifelse(p < 1/2, sqrt(2 * p) - 1, 1 - sqrt(2 * (1-p)))
           },
           epanechnikov = {
             # solve using `polyroot'
             yy <- numeric(n <- length(p))
             yy[p == 0] <- -1
             yy[p == 1] <-  1
             inside <- (p != 0) & (p != 1)
             # coefficients of polynomial (2 + 3 y - y^3)/4
             z <- c(2, 3, 0, -1)/4
             for(i in seq(n)[inside]) {
               sol <- polyroot(z - c(p[i], 0, 0, 0))
               ok <- abs(Im(sol)) < 1e-6
               realpart <- Re(sol)
               ok <- ok & (abs(realpart) <= 1)
               if(sum(ok) != 1)
                 stop(paste("Internal error:", sum(ok), "roots of polynomial"))
               yy[i] <- realpart[ok]
             }
             yy
           },
           biweight = {
             # solve using `polyroot'
             yy <- numeric(n <- length(p))
             yy[p == 0] <- -1
             yy[p == 1] <-  1
             inside <- (p != 0) & (p != 1)
             # coefficients of polynomial (8 + 15 * y - 10 * y^3 + 3 * y^5)/16
             z <- c(8, 15, 0, -10, 0, 3)/16
             for(i in seq(n)[inside]) {
               sol <- polyroot(z - c(p[i], 0, 0, 0, 0, 0))
               ok <- abs(Im(sol)) < 1e-6
               realpart <- Re(sol)
               ok <- ok & (abs(realpart) <= 1)
               if(sum(ok) != 1) 
                 stop(paste("Internal error:", sum(ok), "roots of polynomial"))
               yy[i] <- realpart[ok]
             }
             yy
           },
           cosine = {
             # solve using `uniroot'
             g <- function(y, pval) { (y + sin(pi * y)/pi + 1)/2 - pval }
             yy <- numeric(n <- length(p))
             yy[p == 0] <- -1
             yy[p == 1] <-  1
             inside <- (p != 0) & (p != 1)
             for(i in seq(n)[inside]) 
               yy[i] <- uniroot(g, c(-1,1), pval=p[i])$root
             yy
           },
           optcosine = {
             (2/pi) * asin(2 * p - 1)
           })
  return(mean + a * y)
}

#' integral of t^m k(t) dt from -Inf to r
#'   where k(t) is the standard kernel with support [-1,1]
#' was:    nukernel(r, m, kernel)

kernel.moment <- local({

  kernel.moment <- function(m, r, kernel="gaussian") {
    ker <- match.kernel(kernel)
    check.1.integer(m)
    #' restrict to support
    if(ker != "gaussian") {
      r <- pmin(r, 1)
      r <- pmax(r, -1)
    }
    if(!(m %in% c(0,1,2)) || (ker %in% c("cosine", "optcosine"))) {
      ## use generic integration
      neginf <- if(ker == "gaussian") -10 else -1
      result <- numeric(length(r))
      for(i in seq_along(r))
        result[i] <- integralvalue(kintegrand,
                                   lower=neginf, upper=r[i],
                                   m=m, ker=ker)
      return(result)
    }
    switch(ker,
           gaussian={
             if(m == 0) return(pnorm(r)) else
             if(m == 1) return(-dnorm(r)) else
             return(pnorm(r) - r * dnorm(r))
           },
           rectangular = {
             if(m == 0) return((r + 1)/2) else
             if(m == 1) return((r^2 - 1)/4) else
             return((r^3 + 1)/6)
           },
           triangular={
             m1 <- m+1
             m2 <- m+2
             const <- ((-1)^m1)/m1 + ((-1)^m2)/m2
             answer <- (r^m1)/m1 + ifelse(r < 0, 1, -1) * (r^m2)/m2 - const
             return(answer)
           },
           epanechnikov = {
             if(m == 0)
               return((2 + 3*r - r^3)/4)
             else if(m == 1)
               return((-3 + 6*r^2 - 3*r^4)/16)
             else
               return(( 2 + 5*r^3  - 3* r^5)/20)
           },
           biweight = {
             if(m == 0)
               return((3*r^5 - 10*r^3 + 15*r + 8)/16)
             else if(m == 1)
               return((5*r^6 - 15*r^4 + 15*r^2 -5)/32)
             else 
               return((15*r^7 - 42*r^5 + 35*r^3 + 8)/112)
           },
           # never reached!
           cosine={stop("Sorry, not yet implemented for cosine kernel")},
           optcosine={stop("Sorry, not yet implemented for optcosine kernel")}
           )
  }

  integralvalue <- function(...) integrate(...)$value
  
  kintegrand <- function(x, m, ker) {
    (x^m) * dkernel(x, ker, mean=0, sd=1/kernel.factor(ker))
  }

  kernel.moment
})

kernel.squint <- function(kernel="gaussian", bw=1) {
  kernel <- match.kernel(kernel)
  check.1.real(bw)
  RK <- switch(kernel,
               gaussian = 1/(2 * sqrt(pi)),
               rectangular = sqrt(3)/6, 
               triangular = sqrt(6)/9,
               epanechnikov = 3/(5 * sqrt(5)), 
               biweight = 5 * sqrt(7)/49,
               cosine = 3/4 * sqrt(1/3 - 2/pi^2),
               optcosine = sqrt(1 - 8/pi^2) * pi^2/16)
  return(RK/bw)
}
