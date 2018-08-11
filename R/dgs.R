#
#
#    dgs.R
#
#    $Revision: 1.12 $	$Date: 2018/03/19 14:41:54 $
#
#    Diggle-Gates-Stibbard process
#
#
# -------------------------------------------------------------------
#	

DiggleGatesStibbard <- local({

  # .......... auxiliary functions ................
  dgsTerms <- function(X, Y, idX, idY, rho) {
    stopifnot(is.numeric(rho))
    # sort in increasing order of x coordinate
    oX <- fave.order(X$x)
    oY <- fave.order(Y$x)
    Xsort <- X[oX]
    Ysort <- Y[oY]
    idXsort <- idX[oX]
    idYsort <- idY[oY]
    nX <- npoints(X)
    nY <- npoints(Y)
    # call C routine
    out <- .C("Ediggatsti",
            nnsource = as.integer(nX),
            xsource  = as.double(Xsort$x),
            ysource  = as.double(Xsort$y),
            idsource = as.integer(idXsort),
            nntarget = as.integer(nY),
            xtarget  = as.double(Ysort$x),
            ytarget  = as.double(Ysort$y),
            idtarget = as.integer(idYsort),
            rrho     = as.double(rho),
            values   = as.double(double(nX)),
            PACKAGE = "spatstat")
    answer <- numeric(nX)
    answer[oX] <- out$values
    return(answer)
  }

  # ...... template object ......................
  BlankDGS <- 
    list(
         name   = "Diggle-Gates-Stibbard process",
         creator = "DiggleGatesStibbard",
         family  = "pairwise.family",  # evaluated later
         pot    = function(d, par) {
           rho <- par$rho
           v <- log(sin((pi/2) * d/rho)^2)
           v[ d > par$rho ] <- 0
           attr(v, "IsOffset") <- TRUE
           v
         },
         par    = list(rho = NULL),  # to be filled in later
         parnames = "interaction range",
         hasInf = TRUE,
         init   = function(self) {
           rho <- self$par$rho
           if(!is.numeric(rho) || length(rho) != 1L || rho <= 0)
             stop("interaction range rho must be a positive number")
         },
         update = NULL,       # default OK
         print = NULL,        # default OK
         interpret =  function(coeffs, self) {
           return(NULL)
         },
         valid = function(coeffs, self) {
           return(TRUE)
         },
         project = function(coeffs, self) {
           return(NULL)
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           rho <- self$par$rho
           return(rho)
         },
         version=NULL, # evaluated later
         # fast evaluation is available for the border correction only
         can.do.fast=function(X,correction,par) {
           return(all(correction %in% c("border", "none")))
         },
         fasteval=function(X,U,EqualPairs,pairpot,potpars,correction,
                        splitInf=FALSE, ...) {
           # fast evaluator for DiggleGatesStibbard interaction
           if(!all(correction %in% c("border", "none")))
             return(NULL)
           if(spatstat.options("fasteval") == "test")
             message("Using fast eval for DiggleGatesStibbard")
           dont.complain.about(splitInf)
           rho <- potpars$rho
           idX <- seq_len(npoints(X))
           idU <- rep.int(-1L, npoints(U))
           idU[EqualPairs[,2L]] <- EqualPairs[,1L]
           v <- dgsTerms(U, X, idU, idX, rho)
           v <- matrix(v, ncol=1L)
           attr(v, "IsOffset") <- TRUE
           return(v)
         },
         Mayer=function(coeffs, self) {
           # second Mayer cluster integral
           rho   <- self$par$rho
           return((pi/2 - 2/pi) * rho^2)
         }
         )
  class(BlankDGS) <- "interact"

  DiggleGatesStibbard <- function(rho) {
    instantiate.interact(BlankDGS, list(rho = rho))
  }

  DiggleGatesStibbard <- intermaker(DiggleGatesStibbard, BlankDGS)

  DiggleGatesStibbard
})
