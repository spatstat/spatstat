#	Iest.R
#
#	I function
#
#	$Revision: 1.16 $	$Date: 2019/10/31 03:01:26 $
#
#
#
Iest <- local({

  Iest <- function(X, ...,
                   eps=NULL, r = NULL, breaks = NULL, correction=NULL) {

    X <- as.ppp(X)
    if(!is.multitype(X))
      stop("Only applicable to multitype point patterns")
    marx <- marks(X, dfok=FALSE)
    ntypes <- length(levels(marx))

    Y <- unmark(split(X))
  
    ## relative proportions 
    ni <- sapply(Y, npoints)
    fi <- ni/sum(ni)

    ## J function of pattern regardless of type
    Jdotdot <- Jest(unmark(X),
                    correction=correction, r=r, eps=eps, breaks=breaks, ...)
    rvals <- Jdotdot$r
  
    ## J function of subpattern of each type i
    Jii <- lapply(Y, Jest, r=rvals, correction=correction, eps=eps, ...)
    nrvals <- lengths(lapply(Jii, getElement, name="r"))
    if(length(unique(nrvals)) != 1 || nrvals[1] != length(rvals))
      stop("Internal error: J function objects have different lengths")

    ## initialise fv object
    alim <- attr(Jdotdot, "alim")
    Z <- fv(data.frame(r=rvals, theo=0),
            "r", substitute(I(r), NULL), "theo",
            . ~ r, alim,
            c("r", "%s[pois](r)"),
            c("distance argument r", "theoretical Poisson %s"),
            fname="I")
  
    ## Estimates of each type
    namii <- unlist(lapply(Jii, names))
    namdd <- names(Jdotdot)
    bothnames <- namii[namii %in% namdd]
  
    if("un" %in% bothnames) {
      Jun <- matrix(extract(Jii, "un"), nrow=ntypes, byrow=TRUE)
      Iun <- apply(fi * Jun, 2, sum) - Jdotdot$un
      Z <- bind.fv(Z, data.frame(un=Iun), "hat(%s)[un](r)",
                   "uncorrected estimate of %s", "un")
    }
    if("rs" %in% bothnames) {
      Jrs <- matrix(extract(Jii, "rs"), nrow=ntypes, byrow=TRUE)
      Irs <- apply(fi * Jrs, 2, sum) - Jdotdot$rs    
      Z <- bind.fv(Z, data.frame(rs=Irs), "hat(%s)[rs](r)",
                   "border corrected estimate of %s", "rs")
    }
    if("han" %in% bothnames) {
      Jhan <- matrix(extract(Jii, "han"), nrow=ntypes, byrow=TRUE)
      Ihan <- apply(fi * Jhan, 2, sum) - Jdotdot$han
      Z <- bind.fv(Z, data.frame(han=Ihan), "hat(%s)[han](r)",
                   "Hanisch-style estimate of %s", "han")
    }
    if("km" %in% bothnames) {
      Jkm <- matrix(extract(Jii, "km"), nrow=ntypes, byrow=TRUE)
      Ikm <- apply(fi * Jkm, 2, sum) - Jdotdot$km
      Z <- bind.fv(Z, data.frame(km=Ikm), "hat(%s)[km](r)",
                   "Kaplan-Meier estimate of %s", "km")
    }
    unitname(Z) <- unitname(X)
    attr(Z, "conserve") <- attr(Jdotdot, "conserve")
    return(Z)
  }

  extract <- function(Zlist, what) sapply(Zlist, "[[", i=what)

  Iest
})


