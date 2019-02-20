#
#
#    saturated.S
#
#    $Revision: 1.10 $	$Date: 2019/02/20 03:34:50 $
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
  if(!isTRUE(all.equal(fop, c("d", "par")))
     && !isTRUE(all.equal(fop, c("d", "tx", "tu", "par"))))
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
         hasInf   = NA,
         init     = NULL,
         update   = function(self, ...){
           do.call(Saturated,
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

Saturated <-
    intermaker(Saturated,
               list(creator="Saturated",
                    name="saturated process with user-defined potential",
                    par=formals(Saturated),
                    parnames=list("the potential",
                        "the name of the interaction")))
