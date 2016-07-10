#'
#'   kernel2d.R
#'
#'  Two-dimensional smoothing kernels
#'
#'  $Revision: 1.4 $ $Date: 2016/07/10 08:41:07 $
#'

.Spatstat.2D.KernelTable <- list(
  gaussian=list(
    d = function(x,y, ...) { dnorm(x) * dnorm(y) },
    a = 1),
  Gaussian=list(
    d = function(x,y, ...) { dnorm(x) * dnorm(y) },
    a = 1),
  epanechnikov=list(
    d = function(x,y, ...) { (2/pi) * pmax(1 - (x^2+y^2), 0) },
    a = sqrt(6)),
  Epanechnikov=list(
    d = function(x,y, ...) { (2/pi) * pmax(1 - (x^2+y^2), 0) },
    a = sqrt(6)),
  disc=list(
    d = function(x,y,...) { (1/pi) * as.numeric(x^2 + y^2 <= 1) },
    a = 2)
)

lookup2Dkernel <- function(kernel) {
  nama <- names(.Spatstat.2D.KernelTable)
  if(!is.character(kernel) || is.na(m <- pmatch(kernel, nama))) 
    stop(paste("kernel must be one of the strings",
               commasep(nama, " or "),
               "or a function(x,y) or a pixel image"),
         call.=FALSE)
  return(.Spatstat.2D.KernelTable[[m]])
}

validate2Dkernel <- function(kernel) {
  if(is.im(kernel) || is.function(kernel)) return(TRUE)
  lookup2Dkernel(kernel)
  return(TRUE)
}

evaluate2Dkernel <- function(kernel, x, y, sigma=NULL, varcov=NULL, ...,
                             scalekernel=is.character(kernel)) {

  if(scalekernel) {
    ## kernel adjustment factor 'a'
    a <- 1
    if(is.character(kernel)) {
      info <- lookup2Dkernel(kernel)
      kerfun <- info$d
      a <- info$a  ## kerfun(a * x, a * y) has bandwidth 1
                   ## kerfun() has bandwidth 1/a
    }
    ## transform coordinates to x',y' such that kerfun(x', y')
    ## yields density k(x,y) at desired bandwidth
    if(is.null(varcov)) {
      aos <- a/sigma
      x <- x * aos
      y <- y * aos
      const <- aos^2
    } else {
      SinvH <- matsqrt(solve(varcov))
      XY <- cbind(x, y) %*% (a * SinvH)
      x <- XY[,1]
      y <- XY[,2]
      const <- det(a * SinvH)
    }
  }

  ## now evaluate kernel
  
  if(is.character(kernel)) {
    result <- const * kerfun(x, y)
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

  

  

  
