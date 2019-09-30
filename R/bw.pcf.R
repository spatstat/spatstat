#'
#' bw.pcf.R
#'
#' $Revision: 1.5 $  $Date: 2019/09/30 07:51:52 $
#'
#' bandwidth selection for pcf
#' with least-squares cross-validation method
#' 
#' Original code by: Rasmus Waagepetersen and Abdollah Jalilian
#'
#' References:
#' Guan, Y. (2007). A composite likelihood cross-validation approach in 
#'   selecting bandwidth for the estimation of the pair correlation function. 
#'   Scandinavian Journal of Statistics, 34(2), 336--346. 
#'   DOI: http://doi.org/10.1111/j.1467-9469.2006.00533.x
#' Guan, Y. (2007). A least-squares cross-validation bandwidth 
#'   selection approach in pair correlation function estimations. 
#'   Statistics & Probability Letters, 77(18), 1722--1729. 
#'   DOI: http://doi.org/10.1016/j.spl.2007.04.016

bw.pcf <- function(X, rmax=NULL, lambda=NULL, divisor="r", 
                   kernel="epanechnikov", nr=10000, bias.correct=TRUE, 
                   cv.method=c("compLik", "leastSQ"), simple=TRUE,
                   srange=NULL, ..., verbose=FALSE, warn=TRUE)
{
  stopifnot(is.ppp(X))
  X <- unmark(X)
  win <- Window(X)
  areaW <- area(win)
  nX <- npoints(X)

  cv.method <- match.arg(cv.method)
  kernel <- match.kernel(kernel)
  
  #' maximum distance lag: rmax
  if (is.null(rmax))
    rmax <- rmax.rule("K", win,  nX/areaW)
  if(is.null(srange))
    srange <- c(0, rmax/4)
  #' number of subintervals for discretization of [0, rmax]: nr
  #' length of subintervals
  discr <- rmax / nr
  #' breaks of subintervals
  rs <- seq(0, rmax, length.out= nr + 1)

  #' closepairs distances: \\ u - v \\
  #' Pre-compute close pair distances for use in 'pcf'
  #'   we need close pairs up to a distance rmax + smax
  #'   where 'smax' is the maximum halfwidth of the support of the kernel
  smax <- srange[2] * (if(kernel == "gaussian") 2 else kernel.factor(kernel))
  cpfull <- closepairs(X, rmax + smax, what="all", twice=TRUE)
  
  #' For cross-validation, restrict close pairs to distance rmax 
  ok <- (cpfull$d <= rmax)
  cp <- lapply(cpfull, "[", i=ok)

  ds <- cp$d
  #' determining closepairs distances are in which subinterval
  idx <- round(ds / discr) + 1L
  idx <- pmin.int(idx, nr+1L)
  
  #' translation edge correction factor: /W|/|W \cap W_{u-v}|
  edgewt <- edge.Trans(dx=cp$dx, dy=cp$dy, W=win, paired=TRUE)
  
  if(homogeneous <- is.null(lambda)) {
    #' homogeneous case
    lambda <- nX/areaW
    lambda2area <- lambda^2 * areaW
    pcfargs <- list(X=X, r=rs,
                    divisor=divisor, kernel=kernel, correction="translate",
                    close=cpfull)
    renorm.factor <- 1
  } else {
    # inhomogeneous case: lambda is assumed to be a numeric vector giving
    # the intensity at the points of the point pattern X
    check.nvector(lambda, nX)
    lambda2area <- lambda[cp$i] * lambda[cp$j] * areaW
    pcfargs <- list(X=X, lambda=lambda, r=rs,
                    divisor=divisor, kernel=kernel, correction="translate",
                    close=cpfull)
    renorm.factor <- (areaW/sum(1/lambda))
  }
  
  stuff <- list(cv.method=cv.method,
                kernel=kernel,
                homogeneous=homogeneous,
                bias.correct=bias.correct,
                simple = simple,
                discr=discr,
                rs=rs,
                cp=cp,
                ds=ds,
                idx=idx,
                edgewt=edgewt,
                pcfargs=pcfargs,
                lambda=lambda,
                lambda2area=lambda2area,
                renorm.factor=renorm.factor,
		show=verbose)
  stuff <- list2env(stuff)

  #' find optimum bandwidth
  z <- optimizeWithTrace(CVforPCF, srange, maximum=TRUE, stuff=stuff)

  #' pack up
  ox <- order(z$x)
  sigma  <- z$x[ox]
  cv     <- z$y[ox]
  criterion <- switch(cv.method,
                      compLik = "composite likelihood cross-validation",
                      leastSQ = "least squares cross-validation")
  result <- bw.optim(cv, sigma, which.max(cv),
                     criterion = criterion,
                     warnextreme=warn, hargnames=c("rmax", "srange"),
                     unitname=unitname(X))
  return(result)
}

CVforPCF <- function(bw, stuff) {
  stuff$bw <- bw
  with(stuff, {
    if(show) splat("bw=", bw)
    #' values of pair correlation at breaks of subintervals
    a <- append(pcfargs, list(bw=bw))
    grs <- if(homogeneous) do.call(pcf.ppp, a) else do.call(pcfinhom, a)
    grs <- grs$trans
    #' bias correction
    if (bias.correct) {
      grs <- grs / pkernel(rs, kernel, 0, bw)
      dcorrec <- pkernel(ds, kernel, 0, bw)
    } else {
      dcorrec <- 1
    }
    #' make sure that the estimated pair correlation at origin is finite
    if (!is.finite(grs[1]))
      grs[1] <- grs[2]
    #' approximate the pair correlation values at closepairs distances
    gds <- grs[idx]
    wt <- edgewt / (2 * pi * ds * lambda2area * dcorrec) * renorm.factor
    #' remove pairs to approximate the cross-validation term: g^{-(u, v)}
    if (simple) {
      gds <- gds - 2 * wt * dkernel(0, kernel, 0, bw)
    } else {
      cpi <- cp$i
      cpj <- cp$j
      for (k in 1:length(ds)) {
        exclude <- (cpi == cpi[k]) | (cpj == cpj[k])
        gds[k] <- gds[k] - 2 * sum(wt[exclude] * 
                                   dkernel(ds[k] - ds[exclude],
                                           kernel, 0, bw))
      }
    }
    #' remove negative and zero values
    gds <- pmax.int(.Machine$double.eps, gds)
    switch(cv.method,
           compLik={
             #' composite likelihood cross-validation
             #' the integral term: 2 \pi \int_{0}^{rmax} \hat g(r) r dr
             normconst <- 2 * pi * sum(grs * rs) * discr
             value <- mean(log(gds)) - log(normconst)
           },
           leastSQ={
             #' least squares cross-validation
             #' the integral term: 2 \pi \int_{0}^{rmax} \hat g^2(r) r dr
             normconst <- 2 * pi * sum(grs^2 * rs) * discr
             value <- 2 * sum(gds * edgewt / (lambda2area)) - normconst
           },
           stop("Unrecognised cross-validation method"))
    if(show) splat("value=", value)
    return(value)
  })
}

