#
#
#    multistrhard.S
#
#    $Revision: 2.39 $	$Date: 2018/03/15 07:37:41 $
#
#    The multitype Strauss/hardcore process
#
#    MultiStraussHard()
#                 create an instance of the multitype Strauss/ harcore
#                 point process
#                 [an object of class 'interact']
#	
# -------------------------------------------------------------------
#	

doMultiStraussHard <- local({
  
  # ........  define potential ......................

  MSHpotential <- function(d, tx, tu, par) {
     # arguments:
     # d[i,j] distance between points X[i] and U[j]
     # tx[i]  type (mark) of point X[i]
     # tu[i]  type (mark) of point U[j]
     #
     # get matrices of interaction radii
     r <- par$iradii
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

     # list all UNORDERED pairs of types to be counted
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
     # list all ORDERED pairs of types to be counted
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
       # apply the relevant interaction distance to each pair of points
       rxu <- r[ tx, tu ]
       str <- (d < rxu)
       str[is.na(str)] <- FALSE
       # and the relevant hard core distance
       hxu <- h[ tx, tu ]
       forbid <- (d < hxu)
       forbid[is.na(forbid)] <- FALSE
       # form the potential
       value <- str
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
     return(z)
   }
  # ............... end of potential function ...................

  # .......... auxiliary functions .................
  
  delMSH <- function(which, types, iradii, hradii, ihc) {
    iradii[which] <- NA
    if(any(!is.na(iradii))) {
      # some gamma interactions left
      # return modified MultiStraussHard with fewer gamma parameters
      return(MultiStraussHard(types, iradii, hradii))
    } else if(any(!ihc)) {
      # no gamma interactions left, but some active hard cores
      return(MultiHard(types, hradii))
    } else return(Poisson())
  }

  # ...........................................................
  
  # Set up basic object except for family and parameters

  BlankMSHobject <- 
    list(
         name     = "Multitype Strauss Hardcore process",
         creator  = "MultiStraussHard",
         family   = "pairwise.family", # evaluated later
         pot      = MSHpotential,
         par      = list(types=NULL, iradii=NULL, hradii=NULL),  # to be added
         parnames = c("possible types",
                      "interaction distances",
                      "hardcore distances"),
         pardesc  = c("vector of possible types",
                      "matrix of interaction distances",
                      "matrix of hardcore distances"),
         hasInf   = TRUE,
         selfstart = function(X, self) {
           types <- self$par$types
           hradii <- self$par$hradii
           if(!is.null(types) && !is.null(hradii)) return(self)
           if(is.null(types)) types <- levels(marks(X))
           if(is.null(hradii)) {
             marx <- marks(X)
             d <- nndist(X, by=marx)
             h <- aggregate(d, by=list(from=marx), min)
             h <- as.matrix(h[, -1, drop=FALSE])
             m <- table(marx)
             mm <- outer(m, m, pmin)
             hradii <- h * mm/(mm+1)
             dimnames(hradii) <- list(types, types)
           }
           MultiStraussHard(types=types,hradii=hradii,iradii=self$par$iradii)
	 },
         init     = function(self) {
           types <- self$par$types
           iradii <- self$par$iradii
           hradii <- self$par$hradii
           # hradii could be NULL
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
             MultiPair.checkmatrix(iradii, nt, sQuote("iradii"))
             if(!is.null(hradii))
               MultiPair.checkmatrix(hradii, nt, sQuote("hradii"))
           }
           ina <- is.na(iradii)
           if(all(ina))
             stop(paste("All entries of", sQuote("iradii"),
                        "are NA"))
           if(!is.null(hradii)) {
             hna <- is.na(hradii)
             both <- !ina & !hna
             if(any(iradii[both] <= hradii[both]))
               stop("iradii must be larger than hradii")
           }
         },
         update = NULL,  # default OK
         print = function(self) {
           types <- self$par$types
           iradii <- self$par$iradii
           hradii <- self$par$hradii
           nt <- nrow(iradii)
           if(waxlyrical('gory')) {
             splat(nt, "types of points")
             if(!is.null(types)) {
               splat("Possible types:")
               print(noquote(types))
             } else splat("Possible types:\t not yet determined")
           }
           splat("Interaction radii:")
           dig <- getOption("digits")
           print(signif(iradii, dig))
           if(!is.null(hradii)) {
             splat("Hardcore radii:")
             print(signif(hradii, dig))
           } else splat("Hardcore radii: not yet determined")
           invisible()
         },
        interpret = function(coeffs, self) {
          # get possible types
          typ <- self$par$types
          ntypes <- length(typ)
          # get matrices of interaction radii
          r <- self$par$iradii
          h <- self$par$hradii
          # list all relevant unordered pairs of types
          uptri <- (row(r) <= col(r)) & !is.na(r)
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
           # interaction radii r[i,j]
           iradii <- self$par$iradii
           # hard core radii r[i,j]
           hradii <- self$par$hradii
           # interaction parameters gamma[i,j]
           gamma <- (self$interpret)(coeffs, self)$param$gammas
           # Check that we managed to estimate all required parameters
           required <- !is.na(iradii)
           if(!all(is.finite(gamma[required])))
             return(FALSE)
           # Check that the model is integrable
           # inactive hard cores ...
           ihc <- (is.na(hradii) | hradii == 0)
           # .. must have gamma <= 1
           return(all(gamma[required & ihc] <= 1))
         },
         project = function(coeffs, self) {
           # types
           types <- self$par$types
           # interaction radii r[i,j]
           iradii <- self$par$iradii
           # hard core radii r[i,j]
           hradii <- self$par$hradii
           # interaction parameters gamma[i,j]
           gamma <- (self$interpret)(coeffs, self)$param$gammas
           # required gamma parameters
           required <- !is.na(iradii)
           # active hard cores
           activehard <- !is.na(hradii) & (hradii > 0)
           ihc <- !activehard
           # problems
           gammavalid <- is.finite(gamma) & (activehard | gamma <= 1)
           naughty    <- required & !gammavalid
           if(!any(naughty))
             return(NULL)
           #
           if(spatstat.options("project.fast")) {
             # remove ALL naughty terms simultaneously
             return(delMSH(naughty, types, iradii, hradii, ihc))
           } else {
             # present a list of candidates
             rn <- row(naughty)
             cn <- col(naughty)
             uptri <- (rn <= cn) 
             upn <- uptri & naughty
             rowidx <- as.vector(rn[upn])
             colidx <- as.vector(cn[upn])
#             matindex <- function(v) { matrix(c(v, rev(v)),
#                                              ncol=2, byrow=TRUE) }
             mats <- lapply(as.data.frame(rbind(rowidx, colidx)), matindex)
             inters <- lapply(mats, delMSH,
                              types=types, iradii=iradii,
                              hradii=hradii, ihc=ihc)
             return(inters)           }
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           r <- self$par$iradii
           h <- self$par$hradii
           ractive <- !is.na(r)
           hactive <- !is.na(h)
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
  class(BlankMSHobject) <- "interact"

  matindex <- function(v) { matrix(c(v, rev(v)), ncol=2, byrow=TRUE) }
  
  # Finally define MultiStraussHard function
  doMultiStraussHard <- function(iradii, hradii=NULL, types=NULL) {
    iradii[iradii == 0] <- NA
    if(!is.null(hradii)) hradii[hradii == 0] <- NA
    out <- instantiate.interact(BlankMSHobject,
                                list(types=types,
                                     iradii = iradii, hradii = hradii))
    if(!is.null(types)) {
      dn <- list(types, types)
      dimnames(out$par$iradii) <- dn
      if(!is.null(out$par$hradii)) dimnames(out$par$hradii) <- dn
    }
    return(out)
  }

  doMultiStraussHard
})


MultiStraussHard <- local({

  MultiStraussHard <- function(iradii, hradii, types=NULL) {
    ## try new syntax
    newcall <- match.call()
    newcall[[1]] <- as.name('doMultiStraussHard')
    out <- try(eval(newcall, parent.frame()), silent=TRUE)
    if(is.interact(out))
      return(out)
    ## try old syntax
    oldcall <- match.call(function(types=NULL, iradii, hradii) {})
    oldcall[[1]] <- as.name('doMultiStraussHard')
    out <- try(eval(oldcall, parent.frame()), silent=TRUE)
    if(is.interact(out))
      return(out)
    ## Syntax is wrong: generate error using new syntax rules
    if(missing(hradii)) hradii <- NULL
    doMultiStraussHard(iradii=iradii, hradii=hradii, types=types)
  }


  BlankMSHobject <- get("BlankMSHobject",
                        envir=environment(doMultiStraussHard))
  
  MultiStraussHard <- intermaker(MultiStraussHard, BlankMSHobject)

  MultiStraussHard
})


  
