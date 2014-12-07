#
#
#    ord.S
#
#    $Revision: 1.6 $	$Date: 2014/12/07 10:43:43 $
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

Ord <- local({

  BlankOrd <- 
  list(
         name     = "Ord process with user-defined potential",
         creator  = "Ord",
         family    = "ord.family",
         pot      = NULL,
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
       version=NULL
  )
  class(BlankOrd) <- "interact"

  Ord <- function(pot, name) {
    out <- instantiate.interact(BlankOrd)
    out$pot <- pot
    if(!missing(name)) out$name <- name
  }

  Ord <- intermaker(Ord, BlankOrd)
})


  
