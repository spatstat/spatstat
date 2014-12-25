##
##    hierstrauss.R
##
##    $Revision: 1.2 $	$Date: 2014/12/24 03:15:34 $
##
##    The hierarchical Strauss process
##
##    HierStrauss()    create an instance of the hierarchical Strauss process
##                 [an object of class 'interact']
##	
## -------------------------------------------------------------------
##	

HierStrauss <- local({

  # ......... define interaction potential

  HSpotential <- function(d, tx, tu, par) {
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

     ## list all ordered pairs of types to be checked
     uptri <- par$archy$relation & !is.na(r)
     mark1 <- (lx[row(r)])[uptri]
     mark2 <- (lx[col(r)])[uptri]
     ## corresponding names
     mark1name <- (lxname[row(r)])[uptri]
     mark2name <- (lxname[col(r)])[uptri]
     vname <- apply(cbind(mark1name,mark2name), 1, paste, collapse="x")
     vname <- paste("mark", vname, sep="")
     npairs <- length(vname)
     ## create logical array for result
     z <- array(FALSE, dim=c(dim(d), npairs),
                dimnames=list(character(0), character(0), vname))
     # go....
     if(length(z) > 0) {
       ## assemble the relevant interaction distance for each pair of points
       rxu <- r[ tx, tu ]
       ## apply relevant threshold to each pair of points
       str <- (d <= rxu)
       ## score
       for(i in 1:npairs) {
         # data points with mark m1
         Xsub <- (tx == mark1[i])
         # quadrature points with mark m2
         Qsub <- (tu == mark2[i])
         # assign
         z[Xsub, Qsub, i] <- str[Xsub, Qsub]
       }
     }
     return(z)
   }
   #### end of 'pot' function ####

  # ........ auxiliary functions ..............
  delHS <- function(which, types, radii, archy) {
    radii[which] <- NA
    if(all(is.na(radii))) return(Poisson())
    return(HierStrauss(types=types, radii=radii, archy=archy))
  }
  
  # Set up basic object except for family and parameters
  BlankHSobject <- 
  list(
       name     = "Hierarchical Strauss process",
       creator  = "HierStrauss",
       family   = "hierpair.family", # evaluated later
       pot      = HSpotential,
       par      = list(types=NULL, radii=NULL, archy=NULL), # filled in later
       parnames = c("possible types",
                    "interaction distances",
                    "hierarchical order"),
       selfstart = function(X, self) {
         if(is.null(self$par$types)) types <- levels(marks(X))
         if(is.null(self$par$archy)) archy <- types
         HierStrauss(types=types,radii=self$par$radii,archy=self$par$archy)
       },
       init = function(self) {
         types <- self$par$types
         if(!is.null(types)) {
           radii <- self$par$radii
           nt <- length(types)
           MultiPair.checkmatrix(radii, nt, sQuote("radii"), asymmok=TRUE)
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
         archy <- self$par$archy
         cat(paste(nrow(radii), "types of points\n"))
         if(!is.null(types) && !is.null(archy)) {
           cat("Possible types and ordering: \n")
           print(archy)
         } else if(!is.null(types)) {
           cat("Possible types: \n")
           print(types)
         } else cat("Possible types:\t not yet determined\n")
         cat("Interaction radii:\n")
         print(hiermat(radii, self$par$archy))
         invisible(NULL)
       },
       interpret = function(coeffs, self) {
         # get possible types
         typ <- self$par$types
         ntypes <- length(typ)
         # get matrix of Strauss interaction radii
         r <- self$par$radii
         # list all unordered pairs of types
         uptri <- self$par$archy$relation & !is.na(r)
         index1 <- (row(r))[uptri]
         index2 <- (col(r))[uptri]
         npairs <- length(index1)
         # extract canonical parameters; shape them into a matrix
         gammas <- matrix(NA, ntypes, ntypes)
         dimnames(gammas) <- list(typ, typ)
         gammas[ cbind(index1, index2) ] <- exp(coeffs)
         #
         return(list(param=list(gammas=gammas),
                     inames="interaction parameters gamma_ij",
                     printable=hiermat(round(gammas, 4), self$par$archy)))
       },
       valid = function(coeffs, self) {
         # interaction parameters gamma[i,j]
         gamma <- (self$interpret)(coeffs, self)$param$gammas
         # interaction radii
         radii <- self$par$radii
         # parameters to estimate
         required <- !is.na(radii) & self$par$archy$relation
         # all required parameters must be finite
         if(!all(is.finite(gamma[required]))) return(FALSE)
         # DIAGONAL interaction parameters must be non-explosive
         d <- diag(rep(TRUE, nrow(radii)))
         return(all(gamma[required & d] <= 1))
       },
       project  = function(coeffs, self) {
         # interaction parameters gamma[i,j]
         gamma <- (self$interpret)(coeffs, self)$param$gammas
         # interaction radii and types
         radii <- self$par$radii
         types <- self$par$types
         # problems?
         uptri <- self$par$archy$relation
         required <- !is.na(radii) & uptri
         okgamma  <- !uptri | (is.finite(gamma) & (gamma <= 1))
         naughty  <- required & !okgamma
         # 
         if(!any(naughty))  
           return(NULL)
         if(spatstat.options("project.fast")) {
           # remove ALL naughty terms simultaneously
           return(delHS(naughty, types, radii, archy))
         } else {
           # present a list of candidates
           rn <- row(naughty)
           cn <- col(naughty)
           ord <- self$par$archy$ordering
           uptri <- (ord[rn] <= ord[cn]) 
           upn <- uptri & naughty
           rowidx <- as.vector(rn[upn])
           colidx <- as.vector(cn[upn])
           mats <- lapply(as.data.frame(rbind(rowidx, colidx)),
                          matrix, ncol=2)
           inters <- lapply(mats, delHS, types=types, radii=radii, archy=archy)
           return(inters)
         }
       },
       irange = function(self, coeffs=NA, epsilon=0, ...) {
         r <- self$par$radii
         active <- !is.na(r) & self$par$archy$relation
         if(any(!is.na(coeffs))) {
           gamma <- (self$interpret)(coeffs, self)$param$gammas
           gamma[is.na(gamma)] <- 1
           active <- active & (abs(log(gamma)) > epsilon)
         }
         if(any(active)) return(max(r[active])) else return(0)
       },
       version=NULL # to be added
       )
  class(BlankHSobject) <- "interact"

  # finally create main function
  HierStrauss <- function(radii, types=NULL, archy=NULL) {
    if(!is.null(types)) {
      if(is.null(archy)) archy <- seq_len(length(types))
      archy <- hierarchicalordering(archy, types)
    } 
    out <- instantiate.interact(BlankHSobject,
                                list(types=types,
                                     radii=radii,
                                     archy=archy))
    if(!is.null(types))
      dimnames(out$par$radii) <- list(types, types)
    return(out)
  }

  HierStrauss <- intermaker(HierStrauss, BlankHSobject)
  
  HierStrauss
})


hierarchicalordering <- function(i, s) {
  s <- as.character(s)
  n <- length(s)
  possible <- if(is.character(i)) s else seq_len(n)
  j <- match(i, possible)
  if(any(uhoh <- is.na(j)))
    stop(paste("Unrecognised",
               ngettext(sum(uhoh), "level", "levels"),
               sQuote(i[uhoh]),
               "amongst possible levels",
               commasep(sQuote(s))))
  if(length(j) < n)
    stop("Ordering is incomplete")
  ord <- order(j)
  m <- matrix(, n, n)
  rel <- matrix(ord[row(m)] <= ord[col(m)], n, n)
  dimnames(rel) <- list(s, s)
  x <- list(indices=j, ordering=ord, labels=s, relation=rel)
  class(x) <- "hierarchicalordering"
  x
}

print.hierarchicalordering <- function(x, ...) {
  cat(paste(x$labels[x$indices], collapse=" > "))
  cat("\n")
}
                     
hiermat <- function (x, h) 
{
  stopifnot(is.matrix(x))
  x[] <- as.character(x)
  if(inherits(h, "hierarchicalordering")) ## allows h to be NULL, etc
    x[!(h$relation)] <- ""
  return(noquote(x))
}
