#'
#'   kernel2d.R
#'
#'  Two-dimensional smoothing kernels
#'
#'  $Revision: 1.6 $ $Date: 2016/07/11 02:49:37 $
#'

.Spatstat.2D.KernelTable <- list(
  #' table entries:
  #'   d = density of standardised kernel
  #'   sd = standard deviation of x coordinate, for standardised kernel
  #'   hw = halfwidth of support of standardised kernel (hw=1 for Gaussian)
  gaussian=list(
    d  = function(x,y, ...) { dnorm(x) * dnorm(y) },
    sd  = 1,
    hw = 8),
  Gaussian=list(
    d  = function(x,y, ...) { dnorm(x) * dnorm(y) },
    sd = 1,
    hw = 8),
  epanechnikov=list(
    d  = function(x,y, ...) { (2/pi) * pmax(1 - (x^2+y^2), 0) },
    sd = 1/sqrt(6),
    hw = 1),
  Epanechnikov=list(
    d  = function(x,y, ...) { (2/pi) * pmax(1 - (x^2+y^2), 0) },
    sd = 1/sqrt(6),
    hw = 1),
  disc=list(
    d  = function(x,y,...) { (1/pi) * as.numeric(x^2 + y^2 <= 1) },
    sd = 1/2,
    hw = 1)
)

validate2Dkernel <- function(kernel, fatal=TRUE) {
  if(is.im(kernel) || is.function(kernel)) return(TRUE)
  nama <- names(.Spatstat.2D.KernelTable)
  if(!is.character(kernel) || is.na(m <- pmatch(kernel, nama))) {
    if(!fatal) return(FALSE)
    stop(paste("kernel must be one of the strings",
               commasep(nama, " or "),
               "or a function(x,y) or a pixel image"),
         call.=FALSE)
  }
  return(TRUE)
}

lookup2Dkernel <- function(kernel) {
  validate2Dkernel(kernel)
  if(!is.character(kernel)) return(NULL)
  return(.Spatstat.2D.KernelTable[[kernel]])
}

evaluate2Dkernel <- function(kernel, x, y, sigma=NULL, varcov=NULL, ...,
                             scalekernel=is.character(kernel)) {

  info <- lookup2Dkernel(kernel)

  if(scalekernel) {
    ## kernel adjustment factor 
    sdK <- if(is.character(kernel)) info$sd else 1
    ## transform coordinates to x',y' such that kerfun(x', y')
    ## yields density k(x,y) at desired bandwidth
    if(is.null(varcov)) {
      rr <- sdK/sigma
      x <- x * rr
      y <- y * rr
      const <- rr^2
    } else {
      SinvH <- matsqrt(solve(varcov))
      rSinvH <- sdK * SinvH
      XY <- cbind(x, y) %*% rSinvH
      x <- XY[,1]
      y <- XY[,2]
      const <- det(rSinvH)
    }
  } 

  ## now evaluate kernel
  
  if(is.character(kernel)) {
    kerfun <- info$d
    result <- kerfun(x, y)
    if(scalekernel)
      result <- const * result
    return(result)
  }

  if(is.function(kernel)) {
    argh <- list(...)
    if(length(argh) > 0)
      argh <- argh[names(argh) %in% names(formals(kernel))]
    result <- do.call(kernel, append(list(x, y), argh))
    if(anyNA(result))
      stop("NA values returned from kernel function")
    if(length(result) != length(x))
      stop("Kernel function returned the wrong number of values")
    if(scalekernel)
      result <- const * result
    return(result)
  }

  if(is.im(kernel)) {
    result <- kernel[list(x=x, y=y)]
    if(anyNA(result))
      stop("Domain of kernel image is not large enough")
    return(result)
    if(scalekernel)
      result <- const * result
  } 

  # never reached
  stop("Unrecognised format for kernel")
}

  

  

  
