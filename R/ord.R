#
#
#    ord.S
#
#    $Revision: 1.9 $	$Date: 2020/01/26 04:32:04 $
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
         hasInf   = NA,
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
    return(out)
  }

  Ord <- intermaker(Ord, BlankOrd)
})


  
