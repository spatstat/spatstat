## bandwidth selection with Least-Squres Cross-Validation method
# By: Rasmus Waagepetersen and Abdollah Jalilian
# reference
# Guan, Y. (2007). A composite likelihood cross-validation approach in 
#   selecting bandwidth for the estimation of the pair correlation function. 
#   Scandinavian Journal of Statistics, 34(2), 336–346. 
#   DOI: http://doi.org/10.1111/j.1467-9469.2006.00533.x
# Guan, Y. (2007). A least-squares cross-validation bandwidth 
#   selection approach in pair correlation function estimations. 
#   Statistics & Probability Letters, 77(18), 1722–1729. 
#   DOI: http://doi.org/10.1016/j.spl.2007.04.016

bw.pcf <- function(X, rmax=NULL, lambda=NULL, divisor="r", 
                   kernel="epanechnikov", nr=10000, bias.correct=TRUE, 
                   cv.method="compLik", simple=TRUE, ...)
{
  win <- Window(X)
  areaW <- area(win)
  
  # maximum distance lag: rmax
  if (is.null(rmax))
    rmax <- rmax.rule("K", win,  npoints(X)/areaW)
  #rmax <- rmax + 0.15/sqrt(npoints(X)/areaW)
  # number of subintervals for discretization of [0, rmax]: nr
  # length of subintervals
  discr <- rmax / nr
  # breaks of subintervals
  rs <- seq(0, rmax, length.out= nr + 1)
  
  cp <- closepairs(X, rmax, what="all", twice=TRUE)
  # closepairs distances: \\ u - v \\
  ds <- cp$d
  # determining closepairs distances are in which subinterval
  idx <- round(ds / discr) + 1
  
  # translation edge correction faqctor: /W|/|W \cap W_{u-v}|
  edgewt <- edge.Trans(dx=cp$dx, dy=cp$dy, W=win, paired=TRUE)
  
  if (is.null(lambda))
  {
    # homogeneous case
    lambda <- npoints(X)/areaW
    lambda2area <- lambda^2 * areaW
    gfun <- function(bw)
    {
      pcf(X, r=rs, bw=bw, divisor=divisor, kernel=kernel,
          correction="translate")$trans
    }
    renorm.factor <- 1
  } else
  {
    # inhomogeneous case: lambda is asssumed to be a numeric vector giving
    # the intenisty at the points of the point pattern X
    lambda2area <- lambda[cp$i] * lambda[cp$j] * areaW
    gfun <- function(bw)
    {
      pcfinhom(X, r=rs, lambda=lambda, bw=bw, divisor=divisor, 
               kernel=kernel, correction="translate")$trans
    }
    renorm.factor <- (areaW/sum(1/lambda))
  }
  
  # kernel function
  # integral of kernel function: int_{-h}^{d} k(r) dr
  switch(kernel, epanechnikov={
    kerfun <- function(d, bw)
    {
      h <- sqrt(5) * bw
      ifelse(abs(d) < h, 0.75 * (1 - (d/h)^2)/h, 0)
    }
    integkerfun <- function(d, bw)
    {
      h <- sqrt(5) * bw
      ifelse(d <= h, 3/4 * (1 + d / h) - 1/4 * (1 + d^3 / h^3), 1)
    }
  }, rectangular={
    kerfun <- function(d, bw) 
    {
      h <- sqrt(3) * bw
      ifelse(abs(d) < h, 0.5/h, 0)
    }
    integkerfun <- function(d, bw)
    {
      h <- sqrt(3) * bw
      ifelse(d <= h, 1/2 * (1 + d/h), 1)
    }
  }, stop("the specified kernel is not implemented yet!"))
  
  cvfun <- function(bw)
  {
    # values of pair correlation at breaks of subintervals 
    grs <- gfun(bw)
    # bias correction
    if (bias.correct)
    {
      grs <- grs / integkerfun(rs, bw)
      dcorrec <- integkerfun(ds, bw)
    } else{
      dcorrec <- 1
    }
    # make sure that the estimated pair correlation at origin is finite
    if (!is.finite(grs[1]))
      grs[1] <- grs[2]
    # approximate the pair correlation values at closepairs distances
    gds <- grs[idx]
    wt <- edgewt / (2 * pi * ds * lambda2area * dcorrec) * renorm.factor
    # remove pairs to approximate the cross-validation term: g^{-(u, v)}
    if (simple)
    {
      gds <- gds - 2 * wt * kerfun(0, bw)
    } else{
      for (k in 1:length(ds))
      {
        exclude <- (cp$i == cp$i[k]) | (cp$j == cp$j[k])
        gds[k] <- gds[k] - 2 * sum(wt[exclude] * 
                                     kerfun(ds[k] - ds[exclude], bw))
      }
    }
    
    switch(cv.method, compLik={  # composite likelihood cross-validation
      # the integral term: 2 \pi \int_{0}^{rmax} \hat g(r) r dr
      normconst <- 2 * pi * sum(grs * rs) * discr
      return(mean(log(gds)) - log(normconst))
    }, leastSQ={  # least squares cross-validation
      # the integral term: 2 \pi \int_{0}^{rmax} \hat g^2(r) r dr
      normconst <- 2 * pi * sum(grs^2 * rs) * discr
      return(2 * sum(gds * edgewt / (lambda2area)) - normconst)
    }, stop("wrong cross-validation method"))
  }
  
  dyt <- optimize(cvfun, lower=0, upper=rmax / 4, maximum=TRUE)
  
  return(dyt$maximum)
}

# example
# bw.pcf(redwood)
