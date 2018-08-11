##
##    hierstrhard.R
##
##    $Revision: 1.5 $	$Date: 2018/03/15 07:37:41 $
##
##    The hierarchical Strauss-hard core process
##
## -------------------------------------------------------------------
##	

HierStraussHard <- local({

  # ......... define interaction potential

  HSHpotential <- function(d, tx, tu, par) {
     # arguments:
     # d[i,j] distance between points X[i] and U[j]
     # tx[i]  type (mark) of point X[i]
     # tu[j]  type (mark) of point U[j]
     #
     # get matrices of interaction radii
     r <- par$iradii
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
       # and the relevant hard core distance
       hxu <- h[ tx, tu ]
       forbid <- (d < hxu)
       forbid[is.na(forbid)] <- FALSE
       # form the potential
       value <- str
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
     return(z)
   }
   #### end of 'pot' function ####

  # ........ auxiliary functions ..............
  delHSH <- function(which, types, iradii, hradii, archy, ihc) {
    iradii[which] <- NA
    if(any(!is.na(iradii))) {
      # some gamma interactions left
      # return modified HierStraussHard with fewer gamma parameters
      return(HierStraussHard(types=types, iradii=iradii, hradii=hradii,
                             archy=archy))
    } else if(any(!ihc)) {
      # ihc = inactive hard cores
      # no gamma interactions left, but some active hard cores
      return(HierHard(types=types, hradii=hradii, archy=archy))
    } else return(Poisson())
  }
  
  # Set up basic object except for family and parameters
  BlankHSHobject <- 
  list(
    name     = "Hierarchical Strauss-hard core process",
    creator  = "HierStraussHard",
    family   = "hierpair.family", # evaluated later
    pot      = HSHpotential,
    par      = list(types=NULL, iradii=NULL, hradii=NULL, archy=NULL), 
    parnames = c("possible types",
                 "interaction distances",
                 "hardcore distances",
                 "hierarchical order"),
    pardesc  = c("vector of possible types",
                  "matrix of interaction distances",
                  "matrix of hardcore distances",
                 "hierarchical order"),
    hasInf = TRUE, 
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
        h <- as.matrix(h[, -1L, drop=FALSE])
        m <- table(marx)
        mm <- outer(m, m, pmin)
        hradii <- h * mm/(mm+1)
        dimnames(hradii) <- list(types, types)
        h[!(archy$relation)] <- NA
      }
      HierStraussHard(types=types,hradii=hradii,
                      iradii=self$par$iradii, archy=archy)
    },
    init = function(self) {
      types <- self$par$types
      iradii <- self$par$iradii
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
        MultiPair.checkmatrix(iradii, nt, sQuote("iradii"), asymmok=TRUE)
        if(!is.null(hradii))
          MultiPair.checkmatrix(hradii, nt, sQuote("hradii"), asymmok=TRUE)
      }
      ina <- is.na(iradii)
      if(all(ina))
        stop(paste("All entries of", sQuote("iradii"), "are NA"))
      if(!is.null(hradii)) {
        hna <- is.na(hradii)
        both <- !ina & !hna
        if(any(iradii[both] <= hradii[both]))
          stop("iradii must be larger than hradii")
      }
    },
    update = NULL, # default OK
    print = function(self) {
         iradii <- self$par$iradii
         hradii <- self$par$hradii
         types <- self$par$types
         archy <- self$par$archy
         if(waxlyrical('gory'))
           splat(nrow(iradii), "types of points")
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
         splat("Interaction radii:")
         dig <- getOption("digits")
         print(hiermat(signif(iradii, dig), archy))
         if(!is.null(hradii)) {
           splat("Hardcore radii:")
           print(hiermat(signif(hradii, dig), archy))
         } else splat("Hardcore radii: not yet determined")
         invisible(NULL)
       },
       interpret = function(coeffs, self) {
         # get possible types
         typ <- self$par$types
         ntypes <- length(typ)
         ## get matrices of interaction radii
         r <- self$par$iradii
         h <- self$par$hradii
         ## list all unordered pairs of types
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
                     printable=hiermat(dround(gammas), self$par$archy)))
       },
       valid = function(coeffs, self) {
         # interaction radii r[i,j]
         iradii <- self$par$iradii
         # hard core radii r[i,j]
         hradii <- self$par$hradii
         # interaction parameters gamma[i,j]
         gamma <- (self$interpret)(coeffs, self)$param$gammas
         # parameters to estimate
         required <- !is.na(iradii) & self$par$archy$relation
         # all required parameters must be finite
         if(!all(is.finite(gamma[required]))) return(FALSE)
         # DIAGONAL interactions must be non-explosive
         d <- diag(rep(TRUE, nrow(iradii)))
         activehard <- !is.na(hradii) & (hradii > 0)
         return(all(gamma[required & d & !activehard] <= 1))
       },
       project  = function(coeffs, self) {
         # interaction parameters gamma[i,j]
         gamma <- (self$interpret)(coeffs, self)$param$gammas
         # interaction radii
         iradii <- self$par$iradii
         # hard core radii r[i,j]
         hradii <- self$par$hradii
         types <- self$par$types
         archy <- self$par$archy
         # active hard cores
         activehard <- !is.na(hradii) & (hradii > 0)
         ihc <- !activehard
         # problems?
         uptri <- archy$relation
         required <- !is.na(iradii) & uptri
         offdiag <- !diag(nrow(iradii))
         gammavalid <- is.finite(gamma) & (activehard | offdiag | (gamma <= 1))
         naughty <- required & !gammavalid
         # 
         if(!any(naughty))  
           return(NULL)
         if(spatstat.options("project.fast")) {
           # remove ALL naughty terms simultaneously
           return(delHSH(naughty, types, iradii, hradii, archy, ihc))
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
           inters <- lapply(mats, delHSH, types=types,
                            iradii=iradii, hradii=hradii,
                            archy=archy, ihc=ihc)
           return(inters)
         }
       },
       irange = function(self, coeffs=NA, epsilon=0, ...) {
         r <- self$par$iradii
         h <- self$par$hradii
         ractive <- !is.na(r) & self$par$archy$relation
         hactive <- !is.na(h) & self$par$archy$relation
         if(any(!is.na(coeffs))) {
           gamma <- (self$interpret)(coeffs, self)$param$gammas
           gamma[is.na(gamma)] <- 1
           ractive <- ractive & (abs(log(gamma)) > epsilon)
         }
         if(!any(c(ractive,hactive)))
           return(0)
         else
           return(max(c(r[ractive],h[hactive])))
       },
       version=NULL # to be added
       )
  class(BlankHSHobject) <- "interact"

  # finally create main function
  HierStraussHard <- function(iradii, hradii=NULL, types=NULL, archy=NULL) {
    if(!is.null(types)) {
      if(is.null(archy)) archy <- seq_len(length(types))
      archy <- hierarchicalordering(archy, types)
    }
    iradii[iradii == 0] <- NA
    out <- instantiate.interact(BlankHSHobject,
                                list(types=types,
                                     iradii=iradii,
                                     hradii=hradii,
                                     archy=archy))
    if(!is.null(types)) {
      dn <- list(types, types)
      dimnames(out$par$iradii) <- dn
      if(!is.null(out$par$hradii)) dimnames(out$par$hradii) <- dn
    }
    return(out)
  }

  HierStraussHard <- intermaker(HierStraussHard, BlankHSHobject)
  
  HierStraussHard
})
