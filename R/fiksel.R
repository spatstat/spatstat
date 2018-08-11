#
#
#    fiksel.R
#
#    $Revision: 1.18 $	$Date: 2018/03/15 07:37:41 $
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
      idcount <- as.integer(table(factor(EqualPairs[,2L], levels=1:nU)))
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
            PACKAGE = "spatstat")
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
         hasInf = TRUE, 
         selfstart = function(X, self) {
           # self starter for Fiksel
           nX <- npoints(X)
           if(nX < 2) {
             # not enough points to make any decisions
             return(self)
           }
           md <- minnndist(X)
           if(!is.na(hc <- self$par$hc)) {
             # value fixed by user or previous invocation
             # check it
             if(md < hc)
               warning(paste("Hard core distance is too large;",
                             "some data points will have zero probability"))
             return(self)
           }
           if(md == 0) 
             warning(paste("Pattern contains duplicated points:",
                           "hard core must be zero"))
           # take hc = minimum interpoint distance * n/(n+1)
           hcX <- md * nX/(nX+1)
           Fiksel(r=self$par$r, hc = hcX, kappa=self$par$kappa)
         },
         init   = function(self) {
           r <- self$par$r
           hc <- self$par$hc
           kappa <- self$par$kappa
           check.1.real(r)
           check.1.real(kappa)
           if(!is.na(hc)) {
             check.1.real(hc)
             stopifnot(hc > 0)
             stopifnot(r > hc)
           } else stopifnot(r > 0)
         },
         update = NULL,       # default OK
         print = NULL,         # default OK
         interpret =  function(coeffs, self) {
           a <- as.numeric(coeffs[1L])
           return(list(param=list(a=a),
                       inames="interaction strength a",
                       printable=signif(a)))
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
           if(anyNA(coeffs))
             return(r)
           a <- coeffs[1L]
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
       fasteval=function(X,U,EqualPairs,pairpot,potpars,correction,
                         splitInf=FALSE, ...) {
         ## fast evaluator for Fiksel interaction
         if(!all(correction %in% c("border", "none")))
           return(NULL)
         if(spatstat.options("fasteval") == "test")
           message("Using fast eval for Fiksel")
         r <- potpars$r
         hc <- potpars$hc
         kappa <- potpars$kappa
         hclose <- (strausscounts(U, X, hc, EqualPairs) != 0)
         fikselbit <- fikselterms(U, X, r, kappa, EqualPairs)
         if(!splitInf) {
           answer <- ifelseAX(hclose, -Inf, fikselbit)
           answer <- matrix(answer, ncol=1L)
         } else {
           answer <- fikselbit
           answer <- matrix(answer, ncol=1L)
           attr(answer, "-Inf") <- hclose
         }
         return(answer)
           
       },
       Mayer=function(coeffs, self) {
         # second Mayer cluster integral
         a <- as.numeric(coeffs[1L])
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

  Fiksel <- function(r, hc=NA, kappa) {
    instantiate.interact(BlankFiksel, list(r = r, hc = hc, kappa=kappa))
  }

  Fiksel <- intermaker(Fiksel, BlankFiksel)
  
  Fiksel
})
