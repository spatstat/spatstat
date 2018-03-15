#
#
#    ordthresh.S
#
#    $Revision: 1.12 $	$Date: 2018/03/15 07:37:41 $
#
#    Ord process with threshold potential
#
#    OrdThresh()  create an instance of the Ord process
#                 [an object of class 'interact']
#                 with threshold potential
#	
#
# -------------------------------------------------------------------
#	

OrdThresh <- local({

  BlankOrdThresh <- 
    list(
      name     = "Ord process with threshold potential",
      creator  = "OrdThresh",
      family    = "ord.family",
      pot      = function(d, par) {
        (d <= par$r)
      },
      par      = list(r = NULL),
      parnames = "threshold distance",
      hasInf   = FALSE,
      init     = function(self) {
        r <- self$par$r
        if(!is.numeric(r) || length(r) != 1 || r <= 0)
          stop("threshold distance r must be a positive number")
      },
      update = NULL,  # default OK
      print = NULL,    # default OK
      interpret =  function(coeffs, self) {
        loggamma <- as.numeric(coeffs[1])
        gamma <- exp(loggamma)
        return(list(param=list(gamma=gamma),
                    inames="interaction parameter gamma",
                    printable=dround(gamma)))
      },
      valid = function(coeffs, self) {
        loggamma <- as.numeric(coeffs[1])
        is.finite(loggamma)
      },
      project = function(coeffs, self) {
        if((self$valid)(coeffs, self)) return(NULL) else return(Poisson())
      },
      irange = function(...) {
        return(Inf)
      },
      version=NULL
      )
  class(BlankOrdThresh) <- "interact"

  OrdThresh <- function(r) { instantiate.interact(BlankOrdThresh, list(r=r)) }

  OrdThresh <- intermaker(OrdThresh, BlankOrdThresh)

  OrdThresh
})

