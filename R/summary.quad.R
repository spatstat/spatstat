#
# summary.quad.R
#
#  summary() method for class "quad"
#
#  $Revision: 1.7 $ $Date: 2014/05/08 10:29:25 $
#
summary.quad <- function(object, ..., checkdup=FALSE) {
  verifyclass(object, "quad")
  s <- list(
       data  = summary.ppp(object$data, checkdup=checkdup),
       dummy = summary.ppp(object$dummy, checkdup=checkdup),
       param = object$param)
  doit <- function(ww) {
    if(length(ww) > 0) 
      return(list(range=range(ww), sum=sum(ww)))
    else
      return(NULL)
  }
  w <- object$w
  Z <- is.data(object)
  s$w <- list(all=doit(w), data=doit(w[Z]), dummy=doit(w[!Z]))
  class(s) <- "summary.quad"
  return(s)
}

print.summary.quad <- local({

  summariseweights <- function(ww, blah, dp=3) {
    cat(paste(blah, ":\n\t", sep=""))
    if(is.null(ww)) {
      cat("(None)\n")
      return()
    }
    splat(paste0("range: ",
              "[",
              paste(signif(ww$range, digits=dp), collapse=", "),
              "]\t",
              "total: ",
              signif(ww$sum, digits=dp)))
  }

  print.summary.quad <- function(x, ..., dp=3) {
    cat("Quadrature scheme = data + dummy + weights\n")
    pa <- x$param
    if(is.null(pa))
      cat("created by an unknown function.\n")
    cat("Data pattern:\n")
    print(x$data, dp=dp)
    
    cat("\nDummy quadrature points:\n")
    ## How they were computed
    if(!is.null(pa)) {
      dumpar <- pa$dummy
      if(is.null(dumpar))
        cat("(provided manually)\n")
      else if(!is.null(nd <- dumpar$nd)) {
        splat(paste0("(",
                     if(dumpar$random) "stratified random points in " else NULL,
                     nd[1], " x ", nd[2], " ",
                     if(!dumpar$quasi) "grid" else
                     paste(" =", prod(nd), "quasirandom points"),
                     ", plus 4 corner points)"))
      } else
      cat("(rule for creating dummy points not understood)")
    }
    ## Description of them
    print(x$dummy, dp=dp)

    cat("\nQuadrature weights:\n")
    ## How they were computed
    if(!is.null(pa)) {
      wpar <- pa$weight
      if(is.null(wpar))
        cat("(values provided manually)\n")
      else if(!is.null(wpar$method)) {
        if(wpar$method=="grid") {
          cat(paste("(counting weights based on",
                    wpar$ntile[1], "x", wpar$ntile[2],
                    "array of rectangular tiles)\n"))
        } else if(wpar$method=="dirichlet") {
          cat(paste("(Dirichlet tile areas, computed",
                    if(wpar$exact) "exactly" else "by pixel approximation",
                    ")\n"))
        } else
        cat("(rule for creating dummy points not understood)\n")
      }
    }
    summariseweights(x$w$all, "All weights", dp)
    summariseweights(x$w$data, "Weights on data points", dp)
    summariseweights(x$w$dummy, "Weights on dummy points", dp)

    return(invisible(NULL))
  }

  print.summary.quad
})

    
print.quad <- function(x, ...) {
  cat("Quadrature scheme\n")
  splat(paste(x$data$n, "data points, ", x$dummy$n, "dummy points"))
  cat(paste("Total weight ", sum(x$w), "\n"))
  return(invisible(NULL))
}
