#
#
#    ordthresh.S
#
#    $Revision: 1.10 $	$Date: 2012/01/17 01:19:48 $
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

OrdThresh <- function(r) {
  out <- 
  list(
         name     = "Ord process with threshold potential",
         creator  = "OrdThresh",
         family    = ord.family,
         pot      = function(d, par) {
                         (d <= par$r)
                    },
         par      = list(r = r),
         parnames = "threshold distance",
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
                       printable=round(gamma,4)))
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
       version=versionstring.spatstat()
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}
