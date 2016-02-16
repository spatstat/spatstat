##
##    hierhard.R
##
##    $Revision: 1.2 $	$Date: 2016/02/16 01:39:12 $
##
##    The hierarchical hard core process
##
## -------------------------------------------------------------------
##	

HierHard <- local({

  # ......... define interaction potential

  HHpotential <- function(d, tx, tu, par) {
     # arguments:
     # d[i,j] distance between points X[i] and U[j]
     # tx[i]  type (mark) of point X[i]
     # tu[j]  type (mark) of point U[j]
     #
     # get matrices of interaction radii
     h <- par$hradii
     #
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

     ## list all ordered pairs of types to be checked
     uptri <- par$archy$relation & !is.na(h)
     mark1 <- (lx[row(h)])[uptri]
     mark2 <- (lx[col(h)])[uptri]
     ## corresponding names
     mark1name <- (lxname[row(h)])[uptri]
     mark2name <- (lxname[col(h)])[uptri]
     vname <- apply(cbind(mark1name,mark2name), 1, paste, collapse="x")
     vname <- paste("mark", vname, sep="")
     npairs <- length(vname)
     ## create logical array for result
     z <- array(FALSE, dim=c(dim(d), npairs),
                dimnames=list(character(0), character(0), vname))
     # go....
     if(length(z) > 0) {
       # apply relevant hard core distance to each pair of points
       hxu <- h[ tx, tu ]
       forbid <- (d < hxu)
       forbid[is.na(forbid)] <- FALSE
       # form the potential
       value <- array(0, dim=dim(d))
       value[forbid] <- -Inf
       ## score
       for(i in 1:npairs) {
         # data points with mark m1
         Xsub <- (tx == mark1[i])
         # quadrature points with mark m2
         Qsub <- (tu == mark2[i])
         # assign
         z[Xsub, Qsub, i] <- value[Xsub, Qsub]
       }
     }
     attr(z, "IsOffset") <- TRUE
     return(z)
   }
   #### end of 'pot' function ####

  # Set up basic object except for family and parameters
  BlankHHobject <- 
  list(
    name     = "Hierarchical hard core process",
    creator  = "HierHard",
    family   = "hierpair.family", # evaluated later
    pot      = HHpotential,
    par      = list(types=NULL, hradii=NULL, archy=NULL), 
    parnames = c("possible types",
                 "hardcore distances",
                 "hierarchical order"),
    pardesc  = c("vector of possible types",
                  "matrix of hardcore distances",
                  "hierarchical order"),
    selfstart = function(X, self) {
      types <- self$par$types
      hradii <- self$par$hradii
      archy <- self$par$archy
      if(!is.null(types) && !is.null(hradii) && !is.null(archy)) return(self)
      if(is.null(types)) types <- levels(marks(X))
      if(is.null(archy)) 
        archy <- seq_len(length(types))
      if(!inherits(archy, "hierarchicalordering"))
        archy <- hierarchicalordering(archy, types)
      if(is.null(hradii)) {
        marx <- marks(X)
        d <- nndist(X, by=marx)
        h <- aggregate(d, by=list(from=marx), min)
        h <- as.matrix(h[, -1, drop=FALSE])
        m <- table(marx)
        mm <- outer(m, m, pmin)
        hradii <- h * mm/(mm+1)
        dimnames(hradii) <- list(types, types)
        h[!(archy$relation)] <- NA
      }
      HierHard(types=types,hradii=hradii,archy=archy)
    },
    init = function(self) {
      types <- self$par$types
      hradii <- self$par$hradii
      ## hradii could be NULL
      if(!is.null(types)) {
        if(!is.null(dim(types)))
          stop(paste("The", sQuote("types"),
                     "argument should be a vector"))
        if(length(types) == 0)
          stop(paste("The", sQuote("types"),"argument should be",
                     "either NULL or a vector of all possible types"))
        if(anyNA(types))
          stop("NA's not allowed in types")
        if(is.factor(types)) {
          types <- levels(types)
        } else {
          types <- levels(factor(types, levels=types))
        }
        nt <- length(types)
        if(!is.null(hradii))
          MultiPair.checkmatrix(hradii, nt, sQuote("hradii"), asymmok=TRUE)
      }
    },
    update = NULL, # default OK
    print = function(self) {
         hradii <- self$par$hradii
         types <- self$par$types
         archy <- self$par$archy
         if(waxlyrical('gory'))
           splat(nrow(hradii), "types of points")
         if(!is.null(types) && !is.null(archy)) {
           if(waxlyrical('space')) {
             splat("Possible types and ordering:")
           } else cat("Hierarchy: ")
           print(archy)
         } else if(!is.null(types)) {
           (if(waxlyrical('space')) splat else cat)("Possible types: ")
           print(types)
         } else if(waxlyrical('gory'))
           splat("Possible types:\t not yet determined")
         if(!is.null(hradii)) {
           splat("Hardcore radii:")
           print(hiermat(dround(hradii), archy))
         } else splat("Hardcore radii: not yet determined")
         invisible(NULL)
       },
       interpret = function(coeffs, self) {
        # there are no regular parameters (woo-hoo!)
         return(NULL)
       },
       valid = function(coeffs, self) {
         return(TRUE)
       },
       project  = function(coeffs, self) {
         return(NULL)
       },
       irange = function(self, coeffs=NA, epsilon=0, ...) {
         h <- self$par$hradii
         return(max(0, h, na.rm=TRUE))
       },
       version=NULL # to be added
       )
  class(BlankHHobject) <- "interact"

  # finally create main function
  HierHard <- function(hradii=NULL, types=NULL, archy=NULL) {
    if(!is.null(types)) {
      if(is.null(archy)) archy <- seq_len(length(types))
      archy <- hierarchicalordering(archy, types)
    }
    out <- instantiate.interact(BlankHHobject,
                                list(types=types,
                                     hradii=hradii,
                                     archy=archy))
    if(!is.null(types) && !is.null(out$par$hradii)) 
      dimnames(out$par$hradii) <- list(types,types)
    return(out)
  }

  HierHard <- intermaker(HierHard, BlankHHobject)
  
  HierHard
})
