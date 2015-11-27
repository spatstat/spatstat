#
# summary.quad.R
#
#  summary() method for class "quad"
#
#  $Revision: 1.8 $ $Date: 2015/11/27 06:59:30 $
#

summary.quad <- local({

  sumriz <- function(ww) {
    if(length(ww) > 0) 
      return(list(range=range(ww), sum=sum(ww)))
    else
      return(NULL)
  }

  summary.quad <- function(object, ..., checkdup=FALSE) {
    verifyclass(object, "quad")
    s <- list(
      data  = summary.ppp(object$data, checkdup=checkdup),
      dummy = summary.ppp(object$dummy, checkdup=checkdup),
      param = object$param)
    w <- object$w
    Z <- is.data(object)
    s$w <- list(all   = sumriz(w),
                data  = sumriz(w[Z]),
                dummy = sumriz(w[!Z]))
    class(s) <- "summary.quad"
    return(s)
  }

  summary.quad
})

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
    splat("Quadrature scheme = data + dummy + weights")
    pa <- x$param
    if(is.null(pa))
      splat("created by an unknown function.")
    splat("\nData pattern:")
    print(x$data, dp=dp)
    
    splat("\nDummy quadrature points:")
    ## How they were computed
    if(!is.null(pa)) {
      dumpar <- pa$dummy
      if(is.null(dumpar))
        splat("(provided manually)", indent=5)
      else if(is.character(dmethod <- dumpar$method))
        splat(dmethod, indent=5)
      else if(!is.null(nd <- dumpar$nd)) {
        splat(paste0("(",
                     if(dumpar$random) "stratified random points in " else NULL,
                     nd[1], " x ", nd[2], " ",
                     if(!dumpar$quasi) "grid" else
                     paste(" =", prod(nd), "quasirandom points"),
                     ", plus 4 corner points)"),
              indent=5)
      } else
      splat("(rule for creating dummy points not understood)", indent=5)
    }
    ## Description of them
    print(x$dummy, dp=dp)

    splat("Quadrature weights:")
    ## How they were computed
    if(!is.null(pa)) {
      wpar <- pa$weight
      if(is.null(wpar))
        splat("(values provided manually)", indent=5)
      else if(is.character(wmethod <- wpar$method)) {
        switch(wmethod,
               grid = {
                 splat("(counting weights based on",
                       wpar$ntile[1], "x", wpar$ntile[2],
                       "array of rectangular tiles)",
                       indent=5)
               },
               dirichlet = {
                 splat("(Dirichlet tile areas, computed",
                       if(wpar$exact) "exactly)" else "by pixel approximation)",
                       indent=5)
               },
               splat(wmethod, indent=5)
               )
      } else splat("(rule for creating dummy points not understood)")
    }
    summariseweights(x$w$all, "All weights", dp)
    summariseweights(x$w$data, "Weights on data points", dp)
    summariseweights(x$w$dummy, "Weights on dummy points", dp)

    return(invisible(NULL))
  }

  print.summary.quad
})

    
print.quad <- function(x, ...) {
  splat("Quadrature scheme")
  splat(x$data$n, "data points, ", x$dummy$n, "dummy points")
  splat("Total weight ", sum(x$w))
  return(invisible(NULL))
}
