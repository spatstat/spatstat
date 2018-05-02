#
#
#    strausshard.S
#
#    $Revision: 2.37 $	$Date: 2018/05/02 09:38:36 $
#
#    The Strauss/hard core process
#
#    StraussHard()     create an instance of the Strauss-hardcore process
#                      [an object of class 'interact']
#	
#
# -------------------------------------------------------------------
#	

StraussHard <- local({

  BlankStraussHard <- 
    list(
         name   = "Strauss - hard core process",
         creator = "StraussHard",
         family  = "pairwise.family",  # evaluated later
         pot    = function(d, par) {
           v <- (d <= par$r)
           v[ d <= par$hc ] <-  (-Inf)
           v
         },
         par    = list(r = NULL, hc = NULL), # filled in later
         parnames = c("interaction distance",
                      "hard core distance"), 
         hasInf = TRUE, 
         selfstart = function(X, self) {
           # self starter for StraussHard
           nX <- npoints(X)
           if(nX < 2) {
             # not enough points to make any decisions
             return(self)
           }
           r <- self$par$r
           md <- minnndist(X)
           if(md == 0) {
             warning(paste("Pattern contains duplicated points:",
                           "hard core must be zero"))
             return(StraussHard(r=r, hc=0))
           }
           if(!is.na(hc <- self$par$hc)) {
             # value fixed by user or previous invocation
             # check it
             if(md < hc)
               warning(paste("Hard core distance is too large;",
                             "some data points will have zero probability"))
             return(self)
           }
           # take hc = minimum interpoint distance * n/(n+1)
           hcX <- md * nX/(nX+1)
           StraussHard(r=r, hc = hcX)
         },
         init   = function(self) {
           r <- self$par$r
           hc <- self$par$hc
           if(length(hc) != 1)
             stop("hard core distance must be a single value")
           if(!is.na(hc)) {
             if(!is.numeric(hc) || hc <= 0)
               stop("hard core distance hc must be a positive number, or NA")
             if(!is.numeric(r) || length(r) != 1 || r <= hc)
               stop("interaction distance r must be a number greater than hc")
           }
         },
         update = NULL,       # default OK
         print = NULL,         # default OK
         interpret =  function(coeffs, self) {
           loggamma <- as.numeric(coeffs[1])
           gamma <- exp(loggamma)
           return(list(param=list(gamma=gamma),
                       inames="interaction parameter gamma",
                       printable=dround(gamma)))
         },
         valid = function(coeffs, self) {
           loggamma <- as.numeric(coeffs[1])
           return(is.finite(loggamma))
         },
         project = function(coeffs, self) {
           loggamma <- as.numeric(coeffs[1])
           if(is.finite(loggamma))
             return(NULL)
           hc <- self$par$hc
           if(hc > 0) return(Hardcore(hc)) else return(Poisson()) 
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           r <- self$par$r
           hc <- self$par$hc
           if(anyNA(coeffs))
             return(r)
           loggamma <- coeffs[1]
           if(abs(loggamma) <= epsilon)
             return(hc)
           else
             return(r)
         },
       version=NULL, # evaluated later
       # fast evaluation is available for the border correction only
       can.do.fast=function(X,correction,par) {
         return(all(correction %in% c("border", "none")))
       },
      fasteval=function(X,U,EqualPairs,pairpot,potpars,correction,
                        splitInf=FALSE, ...) {
        #' fast evaluator for StraussHard interaction
        if(!all(correction %in% c("border", "none")))
          return(NULL)
        if(spatstat.options("fasteval") == "test")
          message("Using fast eval for StraussHard")
        r <- potpars$r
        hc <- potpars$hc
        hclose <- (strausscounts(U, X, hc, EqualPairs) != 0)
        rclose <- strausscounts(U, X, r,  EqualPairs)
        if(!splitInf) {
          answer <- ifelseAX(hclose, -Inf, rclose)
          answer <- matrix(answer, ncol=1)
        } else {
          answer <- ifelseAX(hclose, 0, rclose)
          answer <- matrix(answer, ncol=1)
          attr(answer, "-Inf") <- hclose
        }
        return(answer)
      },
       Mayer=function(coeffs, self) {
         # second Mayer cluster integral
         gamma <- exp(as.numeric(coeffs[1]))
         r <- self$par$r
         hc <- self$par$hc
         return(pi * (hc^2 + (1-gamma) * (r^2 - hc^2)))
       },
       delta2 = function(X, inte, correction, ..., sparseOK=FALSE) {
         r  <- inte$par$r
         hc <- inte$par$hc
         #' positive part
         U <- as.ppp(X)
         nU <- npoints(U)
         cl <- weightedclosepairs(U, r, correction=correction, what="indices")

         if(is.null(cl)) # can't handle edge correction
           return(NULL)
         
         v <- sparseMatrix(i=cl$i, j=cl$j, x=cl$weight,
                           dims=c(nU, nU))
         
         #' hard core part
         hcl <- closepairs(U, hc, what="indices")
         ihit <- hcl$i
         jhit <- hcl$j
         vh <- NULL

         if(is.ppp(X)) {
           #' count conflicts between data points
           nhit <- as.integer(table(factor(jhit, levels=seq_len(nU))))
           #' for a conflicting pair X[i], X[j],
           #' status of X[j] will change when X[i] is deleted
           #' iff X[j] is only in conflict with X[i]
           changes <- (nhit == 1)
           if(any(changes)) {
             changesJ <- changes[jhit]
             vh <- sparseMatrix(i=ihit[changesJ], j=jhit[changesJ], x=TRUE,
                                dims=c(nU, nU))
           }
         } else if(is.quad(X)) {
           #' count conflicts with existing data points
           izdat <- is.data(X)
           hitdata <- izdat[ihit]
           nhitdata <- as.integer(table(factor(jhit[hitdata],
                                               levels=seq_len(nU))))
           #' for a conflicting pair U[i], U[j],
           #' status of U[j] will change when U[i] is added/deleted
           #' iff EITHER
           #'     U[i] = X[i] is a data point and
           #'     U[j] is only in conflict with X[i],
           #' OR
           #'     U[i] is a dummy point,
           #'     U[j] has no conflicts with X.
           changesJ <- (hitdata & (nhitdata[jhit] == 1)) |
                       (!hitdata & (nhitdata[jhit] == 0))
           if(any(changesJ)) 
             vh <- sparseMatrix(i=ihit[changesJ], j=jhit[changesJ], x=TRUE,
                                dims=c(nU, nU))
         } else stop("X should be a ppp or quad object")

         # pack up
         if(!sparseOK) {
           v <- as.matrix(v)
           if(!is.null(vh)) vh <- as.matrix(vh)
         }
         attr(v, "deltaInf") <- vh
         return(v)
       }
    )
  class(BlankStraussHard) <- "interact"
  
  StraussHard <- function(r, hc=NA) {
    instantiate.interact(BlankStraussHard, list(r=r, hc=hc))
  }

  StraussHard <- intermaker(StraussHard, BlankStraussHard)
  
  StraussHard
})
