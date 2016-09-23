#
# summary.quad.R
#
#  summary() method for class "quad"
#
#  $Revision: 1.11 $ $Date: 2016/09/23 07:38:07 $
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
    X <- object$data
    D <- object$dummy
    s <- list(
      data  = summary.ppp(X, checkdup=checkdup),
      dummy = summary.ppp(D, checkdup=checkdup),
      param = object$param)
    ## make description of dummy point arrangement
    dpar <- object$param$dummy
    eps.given <- dpar$orig$eps # could be NULL
    eps.actual <- NULL
    if(is.null(dpar)) {
      descrip <- "(provided manually)"
    } else if(is.character(dmethod <- dpar$method)) {
      descrip <- dmethod
    } else if(identical(dpar$quasi, TRUE)) {
      descrip <- paste(npoints(D), "quasirandom dummy points",
                       "plus 4 corner points")
      eps.actual <- 1/(2 * sqrt(intensity(D)))
    } else if(!is.null(nd <- dpar$nd)) {
      nd <- ensure2vector(nd)
      eps.actual <- unique(sidelengths(Frame(D))/nd)
      if(identical(dpar$random, TRUE)) {
        descrip <- paste("systematic random dummy points in",
                         nd[1], "x", nd[2], "grid",
                         "plus 4 corner points")
      } else {
        descrip <- paste(nd[1], "x", nd[2],
                         "grid of dummy points, plus 4 corner points")
      }
    } else descrip <- "(rule for creating dummy points not understood)"
    
    if(!is.null(eps.actual)) {
      uD <- unitname(D)
      s$resolution <- numberwithunit(eps.actual, uD)
      if(!is.null(eps.given)) {
        descrip2 <- paste("dummy spacing:",
                          format(eps.given %unit% uD), "requested,", 
                          format(eps.actual %unit% uD), "actual")
      } else {
        descrip2 <- paste("dummy spacing:", format(eps.actual %unit% uD))
      }
      descrip <- c(descrip, descrip2)
    }
    s$descrip <- descrip
    
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

    parbreak()

    splat("Data pattern:")
    print(x$data, dp=dp)

    parbreak()

    splat("Dummy quadrature points:")
    ## How they were computed
    splat(x$descrip, indent=5)
    parbreak()
    ## What arguments were given
    if(!is.null(orig <- pa$dummy$orig))
      splat("Original dummy parameters:",
            paste0(names(orig), "=", orig, collapse=", "))
    ## Description of the dummy points
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
    if(waxlyrical('extras')) {
      summariseweights(x$w$all, "All weights", dp)
      summariseweights(x$w$data, "Weights on data points", dp)
      summariseweights(x$w$dummy, "Weights on dummy points", dp)
    }
    return(invisible(NULL))
  }

  print.summary.quad
})

print.quad <- function(x, ...) {
  splat("Quadrature scheme")
  splat(x$data$n, "data points,", x$dummy$n, "dummy points")
  if(waxlyrical('extras')) {
    sx <- summary(x)
    splat(sx$descrip, indent=5)
  }
  splat("Total weight", sum(x$w), indent=5)
  return(invisible(NULL))
}


