#
#
#    penttinen.R
#
#    $Revision: 1.1 $	$Date: 2015/08/30 08:21:55 $
#
#    Penttinen pairwise interaction
#
#
# -------------------------------------------------------------------
#	

Penttinen <- local({

  # create blank template object without family and pars

  BlankAntti <-
  list(
       name     = "Penttinen process",
       creator  = "Penttinen",
       family    = "pairwise.family", # evaluated later
       pot      = function(d, par) {
         ans <- numeric(length(d))
         dim(ans) <- dim(d)
         zz <- d/(2 * par$r)
         ok <- (zz < 1)
         z <- zz[ok]
         ans[ok] <- (2/pi) * (acos(z) - z * sqrt(1-z^2))
         return(ans)
       },
       par      = list(r = NULL), # to be filled in
       parnames = "circle radius",
       init     = function(self) {
         r <- self$par$r
         if(!is.numeric(r) || length(r) != 1 || r <= 0)
           stop("interaction distance r must be a positive number")
       },
       update = NULL,  # default OK
       print = NULL,    # default OK
       interpret =  function(coeffs, self) {
         theta <- as.numeric(coeffs[1])
         gamma <- exp(theta)
         return(list(param=list(gamma=gamma),
                     inames="interaction parameter gamma",
                     printable=dround(gamma)))
       },
       valid = function(coeffs, self) {
         theta <- as.numeric(coeffs[1])
         return(is.finite(theta) && (theta <= 0))
       },
       project = function(coeffs, self) {
         if((self$valid)(coeffs, self)) return(NULL) else return(Poisson())
       },
       irange = function(self, coeffs=NA, epsilon=0, ...) {
         r <- self$par$r
         if(any(is.na(coeffs)))
           return(2 * r)
         theta <- coeffs[1]
         if(abs(theta) <= epsilon)
           return(0)
         else
           return(2 * r)
       },
       version=NULL # to be filled in 
       )
  class(BlankAntti) <- "interact"


  # Finally define main function
  
  Penttinen <- function(r) {
    instantiate.interact(BlankAntti, list(r=r))
  }

  Penttinen <- intermaker(Penttinen, BlankAntti)
  
  Penttinen
})

