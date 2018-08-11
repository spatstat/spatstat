#
#
#    lennard.R
#
#    $Revision: 1.22 $	$Date: 2018/03/15 07:37:41 $
#
#    Lennard-Jones potential
#
#
# -------------------------------------------------------------------
#	

LennardJones <- local({

  BlankLJ <- 
    list(
         name     = "Lennard-Jones process",
         creator  = "LennardJones",
         family   = "pairwise.family",  # evaluated later
         pot      = function(d, par) {
           sig0 <- par$sigma0
           if(is.na(sig0)) {
             d6 <- d^{-6}
             p <- array(c(-d6^2,d6),dim=c(dim(d),2))
           } else {
             # expand around sig0 and set large numbers to Inf
             drat <- d/sig0
             d6 <- drat^{-6}
             p <- array(c(-d6^2,d6),dim=c(dim(d),2))
             small <- (drat < 1/4)
             small <- array(c(small, small), dim=c(dim(d), 2))
             p[small] <- -Inf
             big <- (drat > 4)
             big <- array(c(big, big), dim=c(dim(d), 2))
             p[big] <- 0
           }
           return(p)
         },
         par      = list(sigma0=NULL),  # filled in later
         parnames = "Initial approximation to sigma",
         hasInf = TRUE,
         selfstart = function(X, self) {
           # self starter for Lennard Jones
           # attempt to set value of 'sigma0'
           if(!is.na(self$par$sigma0)) {
             # value fixed by user or previous invocation
             return(self)
           }
           if(npoints(X) < 2) {
             # not enough points
             return(self)
           }
           s0 <- minnndist(X)
           if(s0 == 0) {
             warning(paste("Pattern contains duplicated points:",
                           "impossible under Lennard-Jones model"))
             s0 <- mean(nndist(X))
             if(s0 == 0)
               return(self)
           }
           LennardJones(s0)           
         },
         init     = function(...){}, # do nothing
         update = NULL, # default OK
         print = NULL,    # default OK
         interpret =  function(coeffs, self) {
           theta1 <- as.numeric(coeffs[1L])
           theta2 <- as.numeric(coeffs[2L])
           sig0 <- self$par$sigma0
           if(is.na(sig0))
             sig0 <- 1
           if(sign(theta1) * sign(theta2) == 1) {
             sigma <- sig0 * (theta1/theta2)^(1/6)
             epsilon <- (theta2^2)/(4 * theta1)
           } else {
             sigma <- NA
             epsilon <- NA
           }
           return(list(param=list(sigma=sigma, epsilon=epsilon),
                       inames="interaction parameters",
                       printable=signif(c(sigma=sigma,epsilon=epsilon))))
         },
         valid = function(coeffs, self) {
           p <- unlist(self$interpret(coeffs, self)$param)
           return(all(is.finite(p) & (p > 0)))
         },
         project = function(coeffs, self) {
           if((self$valid)(coeffs, self)) return(NULL) else return(Poisson())
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           if(anyNA(coeffs) || epsilon == 0)
             return(Inf)
           sig0 <- self$par$sigma0
           if(is.na(sig0)) sig0 <- 1
           theta1 <- abs(coeffs[1L])
           theta2 <- abs(coeffs[2L])
           return(sig0 * max((theta1/epsilon)^(1/12), (theta2/epsilon)^(1/6)))
         },
       version=NULL # filled in later
  )
  class(BlankLJ) <- "interact"

  LennardJones <- function(sigma0=NA) {
    if(is.null(sigma0) || !is.finite(sigma0))
      sigma0 <- NA
    instantiate.interact(BlankLJ, list(sigma0=sigma0))
  }

  LennardJones <- intermaker(LennardJones, BlankLJ)
  
  LennardJones
})

