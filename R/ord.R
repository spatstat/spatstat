#
#
#    ord.S
#
#    $Revision: 1.5 $	$Date: 2014/10/24 00:22:30 $
#
#    Ord process with user-supplied potential
#
#    Ord()  create an instance of the Ord process
#                 [an object of class 'interact']
#                 with user-supplied potential
#	
#
# -------------------------------------------------------------------
#	

Ord <- function(pot, name) {
  if(missing(name))
    name <- "Ord process with user-defined potential"
  
  out <- 
  list(
         name     = name,
         creator  = "Ord",
         family    = ord.family,
         pot      = pot,
         par      = NULL,
         parnames = NULL,
         init     = NULL,
         update   = function(self, ...){
           do.call(Ord,
                   resolve.defaults(list(...),
                                    list(pot=self$pot, name=self$name)))
         } , 
         print = function(self) {
           cat("Potential function:\n")
           print(self$pot)
           invisible()
         },
       version=versionstring.spatstat()
  )
  class(out) <- "interact"
  return(out)
}
