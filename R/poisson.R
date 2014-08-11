#
#
#    poisson.S
#
#    $Revision: 1.7 $	$Date: 2012/01/16 08:26:08 $
#
#    The Poisson process
#
#    Poisson()    create an object of class 'interact' describing
#                 the (null) interpoint interaction structure
#                 of the Poisson process.
#	
#
# -------------------------------------------------------------------
#	

Poisson <- function() {
  out <- 
  list(
         name     = "Poisson process",
         creator  = "Poisson",
         family   = NULL,
         pot      = NULL,
         par      = NULL,
         parnames = NULL,
         init     = function(...) { },
         update   = function(...) { },
         print    = function(self) {
           cat("Poisson process\n")
           invisible()
         },
         valid = function(...) { TRUE },
         project = function(...) NULL, 
         irange = function(...) { 0 },
         version=versionstring.spatstat()
  )
  class(out) <- "interact"
  return(out)
}
