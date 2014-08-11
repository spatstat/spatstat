#
#
#    saturated.S
#
#    $Revision: 1.6 $	$Date: 2007/01/11 03:36:02 $
#
#    Saturated pairwise process with user-supplied potential
#
#    Saturated()  create a saturated pairwise process
#                 [an object of class 'interact']
#                 with user-supplied potential
#	
#
# -------------------------------------------------------------------
#	

Saturated <- function(pot, name) {
  if(missing(name))
    name <- "Saturated process with user-defined potential"
  
  fop <- names(formals(pot))
  if(!identical(all.equal(fop, c("d", "par")), TRUE)
     && !identical(all.equal(fop, c("d", "tx", "tu", "par")), TRUE))
    stop(paste("Formal arguments of pair potential function",
               sQuote("pot"),
               "must be either (d, par) or (d, tx, tu, par)"))

  out <- 
  list(
         name     = name,
         creator  = "Saturated",
         family    = pairsat.family,
         pot      = pot,
         par      = NULL,
         parnames = NULL,
         init     = NULL,
         update   = function(self, ...){
           do.call(Saturated,
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
