#
#
#    strausshard.S
#
#    $Revision: 2.19 $	$Date: 2013/07/19 02:53:00 $
#
#    The Strauss/hard core process
#
#    StraussHard()     create an instance of the Strauss-hardcore process
#                      [an object of class 'interact']
#	
#
# -------------------------------------------------------------------
#	

StraussHard <- local({

  BlankStraussHard <- 
    list(
         name   = "Strauss - hard core process",
         creator = "StraussHard",
         family  = "pairwise.family",  # evaluated later
         pot    = function(d, par) {
           v <- 1 * (d <= par$r)
           v[ d <= par$hc ] <-  (-Inf)
           v
         },
         par    = list(r = NULL, hc = NULL), # filled in later
         parnames = c("interaction distance",
                      "hard core distance"), 
         selfstart = function(X, self) {
           # self starter for StraussHard
           nX <- npoints(X)
           if(nX < 2) {
             # not enough points to make any decisions
             return(self)
           }
           r <- self$par$r
           md <- min(nndist(X))
           if(md == 0) {
             warning(paste("Pattern contains duplicated points:",
                           "hard core must be zero"))
             return(StraussHard(r=r, hc=0))
           }
           if(!is.na(hc <- self$par$hc)) {
             # value fixed by user or previous invocation
             # check it
             if(md < hc)
               warning(paste("Hard core distance is too large;",
                             "some data points will have zero probability"))
             return(self)
           }
           # take hc = minimum interpoint distance * n/(n+1)
           hcX <- md * nX/(nX+1)
           StraussHard(r=r, hc = hcX)
         },
         init   = function(self) {
           r <- self$par$r
           hc <- self$par$hc
           if(length(hc) != 1)
             stop("hard core distance must be a single value")
           if(!is.na(hc)) {
             if(!is.numeric(hc) || hc <= 0)
               stop("hard core distance hc must be a positive number, or NA")
             if(!is.numeric(r) || length(r) != 1 || r <= hc)
               stop("interaction distance r must be a number greater than hc")
           }
         },
         update = NULL,       # default OK
         print = NULL,         # default OK
         interpret =  function(coeffs, self) {
           loggamma <- as.numeric(coeffs[1])
           gamma <- exp(loggamma)
           return(list(param=list(gamma=gamma),
                       inames="interaction parameter gamma",
                       printable=round(gamma,4)))
         },
         valid = function(coeffs, self) {
           loggamma <- as.numeric(coeffs[1])
           return(is.finite(loggamma))
         },
         project = function(coeffs, self) {
           loggamma <- as.numeric(coeffs[1])
           if(is.finite(loggamma))
             return(NULL)
           hc <- self$par$hc
           if(hc > 0) return(Hardcore(hc)) else return(Poisson()) 
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           r <- self$par$r
           hc <- self$par$hc
           if(any(is.na(coeffs)))
             return(r)
           loggamma <- coeffs[1]
           if(abs(loggamma) <= epsilon)
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
         # fast evaluator for StraussHard interaction
         if(!all(correction %in% c("border", "none")))
           return(NULL)
         if(spatstat.options("fasteval") == "test")
           message("Using fast eval for StraussHard")
         r <- potpars$r
         hc <- potpars$hc
         hclose <- strausscounts(U, X, hc, EqualPairs)
         rclose <- strausscounts(U, X, r,  EqualPairs)
         answer <- ifelseXB(hclose == 0, rclose, -Inf)
         return(matrix(answer, ncol=1))
       },
       Mayer=function(coeffs, self) {
         # second Mayer cluster integral
         gamma <- exp(as.numeric(coeffs[1]))
         r <- self$par$r
         hc <- self$par$hc
         return(pi * (hc^2 + (1-gamma) * (r^2 - hc^2)))
       }
         )
  class(BlankStraussHard) <- "interact"
  
  StraussHard <- function(r, hc=NA) {
    instantiate.interact(BlankStraussHard, list(r=r, hc=hc))
  }

  StraussHard
})
