#'
#'   kernel2d.R
#'
#'  Two-dimensional smoothing kernels
#'
#'  $Revision: 1.14 $ $Date: 2018/09/01 09:09:34 $
#'

.Spatstat.2D.KernelTable <- list(
  #' table entries:
  #'   d = density of standardised kernel
  #'   sd = standard deviation of x coordinate, for standardised kernel
  #'   hw = halfwidth of support of standardised kernel 
  gaussian=list(
    d  = function(x,y, ...) { dnorm(x) * dnorm(y) },
    sd  = 1,
    hw = 8,
    symmetric = TRUE),
  epanechnikov=list(
    d  = function(x,y, ...) { (2/pi) * pmax(1 - (x^2+y^2), 0) },
    sd = 1/sqrt(6),
    hw = 1,
    symmetric = TRUE),
  quartic=list(
    d  = function(x,y, ...) { (3/pi) * pmax(1 - (x^2+y^2), 0)^2 },
    sd = 1/sqrt(8),
    hw = 1,
    symmetric = TRUE),
  disc=list(
    d  = function(x,y,...) { (1/pi) * as.numeric(x^2 + y^2 <= 1) },
    sd = 1/2,
    hw = 1,
    symmetric = TRUE)
)

validate2Dkernel <- function(kernel, fatal=TRUE) {
  if(is.character(match2DkernelName(kernel))) return(TRUE)
  if(is.im(kernel) || is.function(kernel)) return(TRUE)
  if(!fatal) return(FALSE)
  if(is.character(kernel))
    stop(paste("Unrecognised choice of kernel", sQuote(kernel),
               paren(paste("options are",
                           commasep(sQuote(names(.Spatstat.2D.KernelTable)))))),
         call.=FALSE)
  stop(paste("kernel should be a character string,",
             "a pixel image, or a function (x,y)"),
       call.=FALSE)
}

match2DkernelName <- function(kernel) {
  if(!is.character(kernel) || length(kernel) != 1) return(NULL)
  nama <- names(.Spatstat.2D.KernelTable)
  m <- pmatch(kernel, nama)
  if(is.na(m)) return(NULL)
  return(nama[m])
}

lookup2DkernelInfo <- function(kernel) {
  validate2Dkernel(kernel)
  kernel <- match2DkernelName(kernel)
  if(is.null(kernel)) return(NULL)
  return(.Spatstat.2D.KernelTable[[kernel]])
}

evaluate2Dkernel <- function(kernel, x, y, sigma=NULL, varcov=NULL, ...,
                             scalekernel=is.character(kernel)) {

  info <- lookup2DkernelInfo(kernel)

  if(scalekernel) {
    ## kernel adjustment factor 
    sdK <- if(is.character(kernel)) info$sd else 1
    ## transform coordinates to x',y' such that kerfun(x', y')
    ## yields density k(x,y) at desired bandwidth
    if(is.null(varcov)) {
      rr <- sdK/sigma
      x <- x * rr
      y <- y * rr
      scalefactor <- rr^2
    } else {
      SinvH <- matrixinvsqrt(varcov)
      rSinvH <- sdK * SinvH
      XY <- cbind(x, y) %*% rSinvH
      x <- XY[,1L]
      y <- XY[,2L]
      scalefactor <- det(rSinvH)
    }
  } 

  ## now evaluate kernel
  
  if(is.character(kernel)) {
    kerfun <- info$d
    result <- kerfun(x, y)
  } else if(is.function(kernel)) {
    argh <- list(...)
    if(length(argh) > 0)
      argh <- argh[names(argh) %in% names(formals(kernel))]
    result <- do.call(kernel, append(list(x, y), argh))
    if(anyNA(result))
      stop("NA values returned from kernel function")
    if(length(result) != length(x))
      stop("Kernel function returned the wrong number of values")
  } else if(is.im(kernel)) {
    result <- kernel[list(x=x, y=y)]
    if(anyNA(result) || length(result) != length(x))
      stop("Domain of kernel image is not large enough")
  } else stop("Unrecognised format for kernel")
  
  if(scalekernel)
    result <- scalefactor * result

  return(result)
}

cutoff2Dkernel <- function(kernel, sigma=NULL, varcov=NULL, ...,
                           scalekernel=is.character(kernel),
                           cutoff=NULL, fatal=FALSE) {
  info <- lookup2DkernelInfo(kernel)

  ## if scalekernel = FALSE, 'cutoff' is an absolute distance
  ## if scalekernel = TRUE,  'cutoff' is expressed in number of s.d.
  
  if(scalekernel) {
    if(is.null(cutoff)) {
      ## template kernel's standard deviation
      sdK <- info$sd %orifnull% 1
      ## template kernel's halfwidth
      hwK <- info$hw %orifnull% 8
      ## cutoff for kernel with sd=1
      cutoff <- hwK/sdK
    } 
    ## required standard deviation 
    if(!is.null(sigma)) {
      sig <- sigma
    } else if(!is.null(varcov)) {
      lam <- eigen(varcov)$values
      sig <- sqrt(max(lam))
    } else stop("Cannot determine standard deviation")
    ## 
    cutoff <- cutoff * sig
  }

  if(fatal && is.null(cutoff))
    stop(paste("The argument", sQuote("cutoff"), "is required",
               "when a non-Gaussian kernel is specified",
               "and scalekernel=FALSE"),
         call.=FALSE)
  
  return(cutoff)
  
}

  

  
