#
#
#    fiksel.R
#
#    $Revision: 1.8 $	$Date: 2012/07/14 06:36:26 $
#
#    Fiksel interaction 
#    
#    ee Stoyan Kendall Mcke 1987 p 161
#
# -------------------------------------------------------------------
#	

Fiksel <- local({

  # ......... auxiliary functions ...........

  fikselterms <- function(U, X, r, kappa, EqualPairs=NULL) {
    answer <- crossfikselterms(U, X, r, kappa)
    nU <- npoints(U)
    # subtract contrinbutions from identical pairs (exp(-0) = 1 for each)
    if(length(EqualPairs) > 0) {
      idcount <- as.integer(table(factor(EqualPairs[,2], levels=1:nU)))
      answer <- answer - idcount
    }
    return(answer)
  }

  crossfikselterms <- function(X, Y, r, kappa) {
    stopifnot(is.numeric(r))
    # sort in increasing order of x coordinate
    oX <- fave.order(X$x)
    oY <- fave.order(Y$x)
    Xsort <- X[oX]
    Ysort <- Y[oY]
    nX <- npoints(X)
    nY <- npoints(Y)
    # call C routine
    out <- .C("Efiksel",
            nnsource = as.integer(nX),
            xsource  = as.double(Xsort$x),
            ysource  = as.double(Xsort$y),
            nntarget = as.integer(nY),
            xtarget  = as.double(Ysort$x),
            ytarget  = as.double(Ysort$y),
            rrmax    = as.double(r),
            kkappa   = as.double(kappa),
            values   = as.double(double(nX)),
            PACKAGE  = "spatstat")
    answer <- integer(nX)
    answer[oX] <- out$values
    return(answer)
  }


  # ........ template object ..............
  
  BlankFiksel <- 
  list(
         name   = "Fiksel process",
         creator = "Fiksel",
         family  = "pairwise.family",  # evaluated later
         pot    = function(d, par) {
           v <- (d <= par$r) * exp( - d * par$kappa)
           v[ d <= par$hc ] <-  (-Inf)
           v
         },
         par    = list(r = NULL, hc = NULL, kappa=NULL),  # filled in later
         parnames = c("interaction distance",
                      "hard core distance",
                      "rate parameter"), 
         init   = function(self) {
           r <- self$par$r
           hc <- self$par$hc
           kappa <- self$par$kappa
           if(!is.numeric(hc) || length(hc) != 1 || hc <= 0)
             stop("hard core distance hc must be a positive number")
           if(!is.numeric(r) || length(r) != 1 || r <= hc)
             stop("interaction distance r must be a number greater than hardcore dstance hc")
           if(!is.numeric(kappa) || length(kappa) != 1)
             stop("rate parameter kappa must be a single number")
         },
         update = NULL,       # default OK
         print = NULL,         # default OK
         interpret =  function(coeffs, self) {
           a <- as.numeric(coeffs[1])
           return(list(param=list(a=a),
                       inames="interaction strength a",
                       printable=round(a,2)))
         },
         valid = function(coeffs, self) {
           a <- (self$interpret)(coeffs, self)$param$a
           return(is.finite(a))
         },
         project = function(coeffs, self) {
           if((self$valid)(coeffs, self))
             return(NULL)
           hc <- self$par$hc
           if(hc > 0) return(Hardcore(hc)) else return(Poisson()) 
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           r <- self$par$r
           hc <- self$par$hc
           if(any(is.na(coeffs)))
             return(r)
           a <- coeffs[1]
           if(abs(a) <= epsilon)
             return(hc)
           else
             return(r)
         },
       version=NULL, # evaluated later
       # fast evaluation is available for the border correction only
       can.do.fast=function(X,correction,par) {
         return(all(correction %in% c("border", "none")))
       },
       fasteval=function(X,U,EqualPairs,pairpot,potpars,correction, ...) {
         # fast evaluator for Fiksel interaction
         if(!all(correction %in% c("border", "none")))
           return(NULL)
         if(spatstat.options("fasteval") == "test")
           message("Using fast eval for Fiksel")
         r <- potpars$r
         hc <- potpars$hc
         kappa <- potpars$kappa
         hclose <- strausscounts(U, X, hc, EqualPairs)
         fikselbit <- fikselterms(U, X, r, kappa, EqualPairs)
         answer <- ifelse(hclose == 0, fikselbit, -Inf)
         return(matrix(answer, ncol=1))
       },
       Mayer=function(coeffs, self) {
         # second Mayer cluster integral
         a <- as.numeric(coeffs[1])
         r     <- self$par$r
         hc    <- self$par$hc
         kappa <- self$par$kappa
         f <- function(x, kappa, a){ 2 * pi * x *
                                       (1 - exp(a * exp(-x * kappa))) }
         hardbit <- integrate(f=f, lower=hc, upper=r,
                              a=a, kappa=kappa)
         mess <- hardbit[["message"]]
         if(!identical(mess, "OK")) {
           warning(mess)
           return(NA)
         }
         return(pi * hc^2 + hardbit$value)
       }
  )
  class(BlankFiksel) <- "interact"

  Fiksel <- function(r, hc, kappa) {
    instantiate.interact(BlankFiksel, list(r = r, hc = hc, kappa=kappa))
  }

  Fiksel
})
