#
#
#    strauss.R
#
#    $Revision: 2.36 $	$Date: 2016/04/25 02:34:40 $
#
#    The Strauss process
#
#    Strauss()    create an instance of the Strauss process
#                 [an object of class 'interact']
#	
#
# -------------------------------------------------------------------
#	

Strauss <- local({

  # create blank template object without family and pars

  BlankStrauss <-
  list(
       name     = "Strauss process",
       creator  = "Strauss",
       family    = "pairwise.family", # evaluated later
       pot      = function(d, par) {
         d <= par$r
       },
       par      = list(r = NULL), # to be filled in
       parnames = "interaction distance",
       init     = function(self) {
         r <- self$par$r
         if(!is.numeric(r) || length(r) != 1 || r <= 0)
           stop("interaction distance r must be a positive number")
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
         return(is.finite(loggamma) && (loggamma <= 0))
       },
       project = function(coeffs, self) {
         if((self$valid)(coeffs, self)) return(NULL) else return(Poisson())
       },
       irange = function(self, coeffs=NA, epsilon=0, ...) {
         r <- self$par$r
         if(anyNA(coeffs))
           return(r)
         loggamma <- coeffs[1]
         if(abs(loggamma) <= epsilon)
           return(0)
         else
           return(r)
       },
       version=NULL, # to be filled in 
       # fast evaluation is available for the border correction only
       can.do.fast=function(X,correction,par) {
         return(all(correction %in% c("border", "none")))
       },
       fasteval=function(X,U,EqualPairs,pairpot,potpars,correction, ...) {
         # fast evaluator for Strauss interaction
         if(!all(correction %in% c("border", "none")))
           return(NULL)
         if(spatstat.options("fasteval") == "test")
           message("Using fast eval for Strauss")
         r <- potpars$r
         answer <- strausscounts(U, X, r, EqualPairs)
         return(matrix(answer, ncol=1))
       },
       Mayer=function(coeffs, self) {
         # second Mayer cluster integral
         gamma <- exp(as.numeric(coeffs[1]))
         r <- self$par$r
         return((1-gamma) * pi * r^2)
       },
       Percy=function(d, coeffs, par, ...) {
         ## term used in Percus-Yevick type approximation
         gamma <- exp(as.numeric(coeffs[1]))
         R <- par$r
         t <- abs(d/(2*R))
         t <- pmin.int(t, 1)
         y <- 2 * R^2 * (pi * (1-gamma)
                         - (1-gamma)^2 * (acos(t) - t * sqrt(1 - t^2)))
         return(y)
       },
       delta2 = function(X, inte, correction, ..., sparseOK=FALSE) {
         r <- inte$par$r
         X <- as.ppp(X) # algorithm is the same for data and dummy points
         nX <- npoints(X)
	 switch(correction,
	        none = ,
		border = {
                  cl <- closepairs(X, r, what="indices")
		  weight <- 1
		},
		isotropic = ,
		Ripley = {
                  cl <- closepairs(X, r, what="ijd")
		  weight <- edge.Ripley(X[cl$i], cl$d)
		},
		translate = {
                  cl <- closepairs(X, r, what="all")
		  weight <- edge.Trans(dx = cl$dx,
		   	    	       dy = cl$dy,
		                       W = Window(X),
				       paired=TRUE)
		},
		return(NULL)
		)
         v <- sparseMatrix(i=cl$i, j=cl$j, x=as.numeric(weight),
                           dims=c(nX, nX))
         if(!sparseOK)
           v <- as.matrix(v)
         return(v)
       }
       )
  class(BlankStrauss) <- "interact"


  # Finally define main function
  
  Strauss <- function(r) {
    instantiate.interact(BlankStrauss, list(r=r))
  }

  Strauss <- intermaker(Strauss, BlankStrauss)
  
  Strauss
})

# generally accessible functions
      
strausscounts <- function(U, X, r, EqualPairs=NULL) {
  answer <- crosspaircounts(U,X,r)
  # subtract counts of identical pairs
  if(length(EqualPairs) > 0) {
    nU <- npoints(U)
    idcount <- as.integer(table(factor(EqualPairs[,2L], levels=1:nU)))
    answer <- answer - idcount
  }
  return(answer)
}

crosspaircounts <- function(X, Y, r) {
  stopifnot(is.ppp(X))
  stopifnot(is.numeric(r) && length(r) == 1)
  stopifnot(is.finite(r))
  stopifnot(r >= 0)
  # sort in increasing order of x coordinate
  oX <- fave.order(X$x)
  oY <- fave.order(Y$x)
  Xsort <- X[oX]
  Ysort <- Y[oY]
  nX <- npoints(X)
  nY <- npoints(Y)
  # call C routine
  out <- .C("Ccrosspaircounts",
            nnsource = as.integer(nX),
            xsource  = as.double(Xsort$x),
            ysource  = as.double(Xsort$y),
            nntarget = as.integer(nY),
            xtarget  = as.double(Ysort$x),
            ytarget  = as.double(Ysort$y),
            rrmax    = as.double(r),
            counts   = as.integer(integer(nX)),
            PACKAGE = "spatstat")
  answer <- integer(nX)
  answer[oX] <- out$counts
  return(answer)
}

closepaircounts <- function(X, r) {
  stopifnot(is.ppp(X))
  stopifnot(is.numeric(r) && length(r) == 1)
  stopifnot(is.finite(r))
  stopifnot(r >= 0)
  # sort in increasing order of x coordinate
  oX <- fave.order(X$x)
  Xsort <- X[oX]
  nX <- npoints(X)
  # call C routine
  out <- .C("Cclosepaircounts",
            nxy    = as.integer(nX),
            x      = as.double(Xsort$x),
            y      = as.double(Xsort$y),
            rmaxi  = as.double(r),
            counts = as.integer(integer(nX)),
            PACKAGE = "spatstat")
  answer <- integer(nX)
  answer[oX] <- out$counts
  return(answer)
}

