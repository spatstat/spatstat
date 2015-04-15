#
#
#    multistrauss.S
#
#    $Revision: 2.23 $	$Date: 2015/03/31 03:57:11 $
#
#    The multitype Strauss process
#
#    MultiStrauss()    create an instance of the multitype Strauss process
#                 [an object of class 'interact']
#	
# -------------------------------------------------------------------
#	

MultiStrauss <- local({

  # ......... define interaction potential

  MSpotential <- function(d, tx, tu, par) {
     # arguments:
     # d[i,j] distance between points X[i] and U[j]
     # tx[i]  type (mark) of point X[i]
     # tu[j]  type (mark) of point U[j]
     #
     # get matrix of interaction radii r[ , ]
     r <- par$radii
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

     # list all UNORDERED pairs of types to be checked
     # (the interaction must be symmetric in type, and scored as such)
     uptri <- (row(r) <= col(r)) & !is.na(r)
     mark1 <- (lx[row(r)])[uptri]
     mark2 <- (lx[col(r)])[uptri]
     # corresponding names
     mark1name <- (lxname[row(r)])[uptri]
     mark2name <- (lxname[col(r)])[uptri]
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
     # create logical array for result
     z <- array(FALSE, dim=c(dim(d), npairs),
                dimnames=list(character(0), character(0), vname))
     # go....
     if(length(z) > 0) {
       # assemble the relevant interaction distance for each pair of points
       rxu <- r[ tx, tu ]
       # apply relevant threshold to each pair of points
       str <- (d <= rxu)
       # assign str[i,j] -> z[i,j,k] where k is relevant interaction code
       for(i in 1:nordpairs) {
         # data points with mark m1
         Xsub <- (tx == mark1o[i])
         # quadrature points with mark m2
         Qsub <- (tu == mark2o[i])
         # assign
         z[Xsub, Qsub, ucode[i]] <- str[Xsub, Qsub]
       }
     }
     return(z)
   }
   #### end of 'pot' function ####

  # ........ auxiliary functions ..............
  delMS <- function(which, types, radii) {
    radii[which] <- NA
    if(all(is.na(radii))) return(Poisson())
    return(MultiStrauss(types, radii))
  }
  
  # Set up basic object except for family and parameters
  BlankMSobject <- 
  list(
       name     = "Multitype Strauss process",
       creator  = "MultiStrauss",
       family   = "pairwise.family", # evaluated later
       pot      = MSpotential,
       par      = list(types=NULL, radii = NULL), # to be filled in later
       parnames = c("possible types", "interaction distances"),
       pardesc  = c("vector of possible types",
                    "matrix of hardcore distances"),
       selfstart = function(X, self) {
         if(!is.null(self$par$types)) return(self)
         types <- levels(marks(X))
         MultiStrauss(types=types,radii=self$par$radii)
       },
       init = function(self) {
         types <- self$par$types
         if(!is.null(types)) {
           radii <- self$par$radii
           nt <- length(types)
           MultiPair.checkmatrix(radii, nt, sQuote("radii"))
           if(length(types) == 0)
             stop(paste("The", sQuote("types"),"argument should be",
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
       update = NULL, # default OK
       print = function(self) {
         radii <- self$par$radii
         types <- self$par$types
         if(waxlyrical('gory')) {
           splat(nrow(radii), "types of points")
           if(!is.null(types)) {
             splat("Possible types: ")
             print(noquote(types))
           } else splat("Possible types:\t not yet determined")
         }
         cat("Interaction radii:\n")
         print(signif(radii, getOption("digits")))
         invisible()
       },
       interpret = function(coeffs, self) {
         # get possible types
         typ <- self$par$types
         ntypes <- length(typ)
         # get matrix of Strauss interaction radii
         r <- self$par$radii
         # list all unordered pairs of types
         uptri <- (row(r) <= col(r)) & (!is.na(r))
         index1 <- (row(r))[uptri]
         index2 <- (col(r))[uptri]
         npairs <- length(index1)
         # extract canonical parameters; shape them into a matrix
         gammas <- matrix(, ntypes, ntypes)
         dimnames(gammas) <- list(typ, typ)
         expcoef <- exp(coeffs)
         gammas[ cbind(index1, index2) ] <- expcoef
         gammas[ cbind(index2, index1) ] <- expcoef
         #
         return(list(param=list(gammas=gammas),
                     inames="interaction parameters gamma_ij",
                     printable=dround(gammas)))
       },
       valid = function(coeffs, self) {
         # interaction parameters gamma[i,j]
         gamma <- (self$interpret)(coeffs, self)$param$gammas
         # interaction radii
         radii <- self$par$radii
         # parameters to estimate
         required <- !is.na(radii)
         gr <- gamma[required]
         return(all(is.finite(gr) & gr <= 1))
       },
       project  = function(coeffs, self) {
         # interaction parameters gamma[i,j]
         gamma <- (self$interpret)(coeffs, self)$param$gammas
         # interaction radii and types
         radii <- self$par$radii
         types <- self$par$types
         # problems?
         required <- !is.na(radii)
         okgamma  <- is.finite(gamma) & (gamma <= 1)
         naughty  <- required & !okgamma
         # 
         if(!any(naughty))  
           return(NULL)
         if(spatstat.options("project.fast")) {
           # remove ALL naughty terms simultaneously
           return(delMS(naughty, types, radii))
         } else {
           # present a list of candidates
           rn <- row(naughty)
           cn <- col(naughty)
           uptri <- (rn <= cn) 
           upn <- uptri & naughty
           rowidx <- as.vector(rn[upn])
           colidx <- as.vector(cn[upn])
           matindex <- function(v) { matrix(c(v, rev(v)),
                                            ncol=2, byrow=TRUE) }
           mats <- lapply(as.data.frame(rbind(rowidx, colidx)), matindex)
           inters <- lapply(mats, delMS, types=types, radii=radii)
           return(inters)
         }
       },
       irange = function(self, coeffs=NA, epsilon=0, ...) {
         r <- self$par$radii
         active <- !is.na(r)
         if(any(!is.na(coeffs))) {
           gamma <- (self$interpret)(coeffs, self)$param$gammas
           gamma[is.na(gamma)] <- 1
           active <- active & (abs(log(gamma)) > epsilon)
         }
         if(any(active)) return(max(r[active])) else return(0)
       },
       version=NULL # to be added
       )
  class(BlankMSobject) <- "interact"

  # finally create main function
  MultiStrauss <- function(radii, types=NULL) {
    if((missing(radii) || !is.matrix(radii)) && is.matrix(types)) {
      ## old syntax: (types=NULL, radii)
      radii <- types
      types <- NULL
    }
    radii[radii == 0] <- NA
    out <- instantiate.interact(BlankMSobject, list(types=types, radii = radii))
    if(!is.null(types))
      dimnames(out$par$radii) <- list(types, types)
    return(out)
  }

  MultiStrauss <- intermaker(MultiStrauss, BlankMSobject)
  
  MultiStrauss
})
