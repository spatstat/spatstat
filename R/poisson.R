#
#
#    poisson.S
#
#    $Revision: 1.8 $	$Date: 2015/10/21 09:06:57 $
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

Poisson <- local({

  BlankPoisson <- list(
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
    version=NULL
    )
  
  class(BlankPoisson) <- "interact"

  Poisson <- function() { BlankPoisson }

  Poisson <- intermaker(Poisson, BlankPoisson)

  Poisson
})
                 
