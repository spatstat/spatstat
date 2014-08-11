#
#
#    hardcore.S
#
#    $Revision: 1.7 $	$Date: 2012/06/28 04:20:25 $
#
#    The Hard core process
#
#    Hardcore()     create an instance of the Hard Core process
#                      [an object of class 'interact']
#	
#
# -------------------------------------------------------------------
#	

Hardcore <- local({

  BlankHardcore <- 
  list(
         name   = "Hard core process",
         creator = "Hardcore",
         family  = "pairwise.family",  # evaluated later
         pot    = function(d, par) {
           v <- 0 * d
           v[ d <= par$hc ] <-  (-Inf)
           attr(v, "IsOffset") <- TRUE
           v
         },
         par    = list(hc = NULL),  # filled in later
         parnames = "hard core distance", 
         init   = function(self) {
           hc <- self$par$hc
           if(!is.numeric(hc) || length(hc) != 1 || hc <= 0)
             stop("hard core distance hc must be a positive number")
         },
         update = NULL,       # default OK
         print = NULL,        # default OK
         interpret =  function(coeffs, self) {
           return(NULL)
         },
         valid = function(coeffs, self) {
           return(TRUE)
         },
         project = function(coeffs, self) {
           return(NULL)
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           hc <- self$par$hc
           return(hc)
         },
       version=NULL, # evaluated later
       # fast evaluation is available for the border correction only
       can.do.fast=function(X,correction,par) {
         return(all(correction %in% c("border", "none")))
       },
       fasteval=function(X,U,EqualPairs,pairpot,potpars,correction, ...) {
         # fast evaluator for Hardcore interaction
         if(!all(correction %in% c("border", "none")))
           return(NULL)
         if(spatstat.options("fasteval") == "test")
           message("Using fast eval for Hardcore")
         hc <- potpars$hc
         # call evaluator for Strauss process
         counts <- strausscounts(U, X, hc, EqualPairs)
         # all counts should be zero
         v <- matrix(ifelse(counts > 0, -Inf, 0), ncol=1)
         attr(v, "IsOffset") <- TRUE
         return(v)
       },
       Mayer=function(coeffs, self) {
         # second Mayer cluster integral
         hc <- self$par$hc
         return(pi * hc^2)
       }
  )
  class(BlankHardcore) <- "interact"
  
  Hardcore <- function(hc) {
    instantiate.interact(BlankHardcore, list(hc=hc))
  }

  Hardcore
})
