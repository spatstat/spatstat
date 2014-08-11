#
#
#    multihard.R
#
#    $Revision: 1.9 $	$Date: 2013/05/01 10:17:27 $
#
#    The Hard core process
#
#    Hardcore()     create an instance of the Hard Core process
#                      [an object of class 'interact']
#	
#
# -------------------------------------------------------------------
#	

MultiHard <- local({

  # .... multitype hard core potential
  
  MHpotential <- function(d, tx, tu, par) {
     # arguments:
     # d[i,j] distance between points X[i] and U[j]
     # tx[i]  type (mark) of point X[i]
     # tu[i]  type (mark) of point U[j]
     #
     # get matrices of interaction radii
     h <- par$hradii

     # get possible marks and validate
     if(!is.factor(tx) || !is.factor(tu))
	stop("marks of data and dummy points must be factor variables")
     lx <- levels(tx)
     lu <- levels(tu)
     if(length(lx) != length(lu) || any(lx != lu))
	stop("marks of data and dummy points do not have same possible levels")

     if(!identical(lx, par$types))
        stop("data and model do not have the same possible levels of marks")
     if(!identical(lu, par$types))
        stop("dummy points and model do not have the same possible levels of marks")

     # ensure factor levels are acceptable for column names (etc)
     lxname <- make.names(lx, unique=TRUE)
     
     # list all UNORDERED pairs of types to be checked
     # (the interaction must be symmetric in type, and scored as such)
     uptri <- (row(h) <= col(h)) & (!is.na(h))
     mark1 <- (lx[row(h)])[uptri]
     mark2 <- (lx[col(h)])[uptri]
     # corresponding names
     mark1name <- (lxname[row(h)])[uptri]
     mark2name <- (lxname[col(h)])[uptri]
     vname <- apply(cbind(mark1name,mark2name), 1, paste, collapse="x")
     vname <- paste("mark", vname, sep="")
     npairs <- length(vname)
     # list all ORDERED pairs of types to be checked
     # (to save writing the same code twice)
     different <- mark1 != mark2
     mark1o <- c(mark1, mark2[different])
     mark2o <- c(mark2, mark1[different])
     nordpairs <- length(mark1o)
     # unordered pair corresponding to each ordered pair
     ucode <- c(1:npairs, (1:npairs)[different])
     #
     # create numeric array for result
     z <- array(0, dim=c(dim(d), npairs),
                dimnames=list(character(0), character(0), vname))
     # go....
     if(length(z) > 0) {
       # apply the relevant hard core distance to each pair of points
       hxu <- h[ tx, tu ]
       forbid <- (d < hxu)
       forbid[is.na(forbid)] <- FALSE
       # form the potential 
       value <- array(0, dim=dim(d))
       value[forbid] <- -Inf
       # assign value[i,j] -> z[i,j,k] where k is relevant interaction code
       for(i in 1:nordpairs) {
         # data points with mark m1
         Xsub <- (tx == mark1o[i])
         # quadrature points with mark m2
         Qsub <- (tu == mark2o[i])
         # assign
         z[Xsub, Qsub, ucode[i]] <- value[Xsub, Qsub]
       }
     }
     attr(z, "IsOffset") <- TRUE
     return(z)
   }
   #### end of 'pot' function ####

  # ............ template object ...................
  
  BlankMH <- 
  list(
       name     = "Multitype Hardcore process",
       creator  = "MultiHard",
       family   = "pairwise.family",  # evaluated later
       pot      = MHpotential,
       par      = list(types=NULL, hradii = NULL), # filled in later
       parnames = c("possible types", "hardcore distances"),
       selfstart = function(X, self) {
         if(!is.null(self$par$types)) return(self)
         types <- levels(marks(X))
         MultiHard(types=types,hradii=self$par$hradii)
       },
       init     = function(self) {
         types <- self$par$types
         if(!is.null(types)) {
           h <- self$par$hradii
           nt <- length(types)
           MultiPair.checkmatrix(h, nt, sQuote("hradii"))
           if(length(types) == 0)
             stop(paste("The", sQuote("types"),
                        "argument should be",
                        "either NULL or a vector of all possible types"))
           if(any(is.na(types)))
             stop("NA's not allowed in types")
           if(is.factor(types)) {
             types <- levels(types)
           } else {
             types <- levels(factor(types, levels=types))
           }
         }
       },
       update = NULL,  # default OK
       print = function(self) {
         print.isf(self$family)
         cat(paste("Interaction:\t", self$name, "\n"))
         h <- self$par$hradii
         cat(paste(nrow(h), "types of points\n"))
         types <- self$par$types
         if(!is.null(types)) {
           cat("Possible types: \n")
           print(types)
         } else cat("Possible types: \t not yet determined\n")
         cat("Hardcore radii:\n")
         print(h)
         invisible()
       },
       interpret = function(coeffs, self) {
        # there are no regular parameters (woo-hoo!)
         return(NULL)
       },
       valid = function(coeffs, self) {
         return(TRUE)
       },
       project = function(coeffs, self) {
         return(NULL)
       },
       irange = function(self, coeffs=NA, epsilon=0, ...) {
         h <- self$par$hradii
         return(max(0, h, na.rm=TRUE))
       },
       version=NULL # fix later
  )
  class(BlankMH) <- "interact"

  MultiHard <- function(types=NULL, hradii) {
    out <- instantiate.interact(BlankMH, list(types=types, hradii = hradii))
    if(!is.null(types))
      dimnames(out$par$hradii) <- list(types, types)
    return(out)
  }

  MultiHard
})
