#
#
#    ord.S
#
#    $Revision: 1.4 $	$Date: 2007/01/11 03:36:02 $
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
           cat(paste(self$name, "\n"))
           cat("Potential function:\n")
           print(self$pot)
           invisible()
         },
       version=versionstring.spatstat()
  )
  class(out) <- "interact"
  return(out)
}
