#
#
#    areainter.R
#
#    $Revision: 1.17 $	$Date: 2012/01/18 09:37:46 $
#
#    The area interaction
#
#    AreaInter()    create an instance of the area-interaction process
#                 [an object of class 'interact']
#	
#
# -------------------------------------------------------------------
#

AreaInter <- local({

  # area-interaction potential function
  areapot <- 
    function(X,U,EqualPairs,pars,correction, ...) {
      if(any(correction != "border"))
        warning(paste("Only the border correction is available;",
                      "other options were ignored"))
      n <- U$n
      answer <- numeric(n)
      r <- pars$r
      if(is.null(r)) stop("internal error: r parameter not found")
      dummies <- !(seq_len(n) %in% EqualPairs[,2])
      if(sum(dummies) > 0)
        answer[dummies] <- -areaGain(U[dummies], X, r)
      ii <- EqualPairs[,1]
      jj <- EqualPairs[,2]
      answer[jj] <- -areaLoss(X, r, subset=ii)
#    for(k in seq_len(nrow(EqualPairs))) {
#      i <- EqualPairs[k,1]
#      j <- EqualPairs[k,2]
#      answer[j] <- -areaGain(U[j], X[-i], r)
#    }
      return(1 + answer/(pi * r^2))
    }

  # template object without family, par, version
  BlankAI <- 
  list(
         name     = "Area-interaction process",
         creator  = "AreaInter",
         family   = "inforder.family", # evaluated later
         pot      = areapot,
         par      = list(r = NULL), # to be filled in
         parnames = "disc radius",
         init     = function(self) {
                      r <- self$par$r
                      if(!is.numeric(r) || length(r) != 1 || r <= 0)
                       stop("disc radius r must be a positive number")
                    },
         update = NULL,  # default OK
         print = NULL,    # default OK
         interpret =  function(coeffs, self) {
           logeta <- as.numeric(coeffs[1])
           eta <- exp(logeta)
           return(list(param=list(eta=eta),
                       inames="interaction parameter eta",
                       printable=round(eta,4)))
         },
         valid = function(coeffs, self) {
           eta <- ((self$interpret)(coeffs, self))$param$eta
           return(is.finite(eta))
         },
         project = function(coeffs, self) {
           if((self$valid)(coeffs, self))
             return(NULL)
           return(Poisson())
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           r <- self$par$r
           if(any(is.na(coeffs)))
             return(2 * r)
           logeta <- coeffs[1]
           if(abs(logeta) <= epsilon)
             return(0)
           else
             return(2 * r)
         },
       version=NULL # to be added
  )
  class(BlankAI) <- "interact"

  AreaInter <- function(r) {
    instantiate.interact(BlankAI, list(r=r))
  }

  AreaInter
})


             
