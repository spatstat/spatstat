#
#
#    hybrid.R
#
#    $Revision: 1.6 $	$Date: 2015/10/21 09:06:57 $
#
#    Hybrid of several interactions
#
#    Hybrid()    create a hybrid of several interactions
#                 [an object of class 'interact']
#	
#
# -------------------------------------------------------------------
#	

Hybrid <- local({

  Hybrid <- function(...) {
    interlist <- list(...)
    n <- length(interlist)
    if(n == 0)
      stop("No arguments given")
    #' arguments may be interaction objects or ppm objects
    isinter <- unlist(lapply(interlist, is.interact))
    isppm   <- unlist(lapply(interlist, is.ppm))
    if(any(nbg <- !(isinter | isppm)))
      stop(paste(ngettext(sum(nbg), "Argument", "Arguments"),
                 paste(which(nbg), collapse=", "),
                 ngettext(sum(nbg), "is not an interaction",
                          "are not interactions")))
    #' ensure the list contains only interaction objects
    if(any(isppm))
      interlist[isppm] <- lapply(interlist[isppm], as.interact)
    #' recursively expand any components that are themselves hybrids
    while(any(ishybrid <- unlist(lapply(interlist, is.hybrid)))) {
      i <- min(which(ishybrid))
      n <- length(interlist)
      expandi <- interlist[[i]]$par
      interlist <- c(if(i > 1) interlist[1:(i-1)] else NULL,
                     expandi,
                     if(i < n) interlist[(i+1):n] else NULL)
    }
    #' 
    ncomponents <- length(interlist)
    if(ncomponents == 1) {
      #' single interaction - return it
      return(interlist[[1]])
    }
    #' ensure all components have names
    names(interlist) <- good.names(names(interlist),
                                   "HybridComponent", 1:ncomponents)
    out <- 
      list(
        name     = "Hybrid interaction",
        creator  = "Hybrid",
        family    = hybrid.family,
        pot      = NULL,
        par      = interlist,
        parnames = names(interlist),
        selfstart = function(X, self) {
          ilist <- self$par
          sslist <- lapply(ilist, getElement, name="selfstart")
          has.ss <- sapply(sslist, is.function)
          if(any(has.ss)) {
            ilist[has.ss] <- lapply(ilist[has.ss], invokeSelfStart, Y=X)
            self$par <- ilist
          }
          return(self)
        },
        init     = NULL,
        update = NULL,  # default OK
        print = function(self, ..., family=FALSE, brief=FALSE) {
          if(family)
            print.isf(self$family)
          ncomponents <- length(self$par)
          clabs <- self$parnames
          splat("Hybrid of", ncomponents, "components:",
                commasep(sQuote(clabs)))
          for(i in 1:ncomponents) {
            splat(paste0(clabs[i], ":"))
            print(self$par[[i]], ..., family=family, brief=brief)
          }
          parbreak()
          return(invisible(NULL))
        },
        interpret =  function(coeffs, self) {
          interlist <- self$par
          result <- list(param=list(),
                         inames=character(0),
                         printable=list())
          for(i in 1:length(interlist)) {
            interI <- interlist[[i]]
            nameI  <- names(interlist)[[i]]
            nameI. <- paste(nameI, ".", sep="")
            #' find coefficients with prefix that exactly matches nameI.
            Cname  <- names(coeffs)
            prefixlength <- nchar(nameI.)
            Cprefix <- substr(Cname, 1, prefixlength)
            relevant <- (Cprefix == nameI.)
            #' extract them
            if(any(relevant)) {
              Crelevant <- coeffs[relevant]
              names(Crelevant) <-
                substr(Cname[relevant], prefixlength+1, max(nchar(Cname)))
              #' invoke the self-interpretation of interI
              interpretI <- interI$interpret
              if(is.function(interpretI)) {
                resultI <- interpretI(Crelevant, interI)
                paramI  <- resultI$param
                prinI   <- resultI$printable
                inamesI <- resultI$inames
                inamesI <- paste(nameI, inamesI)
                if(length(prinI) > 0) {
                  result$param     <- append(result$param, paramI)
                  result$printable <- append(result$printable, list(prinI))
                  result$inames <- c(result$inames, inamesI)
                }
              }
            }
          }
          return(result)
        },
        valid = function(coeffs, self) {
          #' check validity via mechanism used for 'rmhmodel' 
          siminfo <- .Spatstat.Rmhinfo[["Hybrid interaction"]]
          Z <- siminfo(coeffs, self)
          cifs   <- Z$cif
          pars   <- Z$par
          ntypes <- Z$ntypes
          if((Ncif <- length(cifs)) == 1) {
            #' single cif
            pars <- append(pars, list(beta=rep.int(1, ntypes)))
          } else {
            for(i in 1:Ncif) 
              pars[[i]] <- append(pars[[i]], list(beta=rep.int(1, ntypes[i])))
          }
          RM <- rmhmodel(cif=cifs, par=pars, types=1:max(ntypes), 
                         stopinvalid=FALSE)
          return(RM$integrable)
        },
        project = function(coeffs, self) {
          if((self$valid)(coeffs, self)) return(NULL)
          #' separate into components
          spl <- splitHybridInteraction(coeffs, self)
          interlist <- spl$interlist
          coeflist  <- spl$coeflist
          #' compute projection for each component interaction
          Ncif <- length(interlist)
          projlist <- vector(mode="list", length=Ncif)
          nproj    <- integer(Ncif)
          for(i in 1:Ncif) {
            coefsI <- coeflist[[i]]
            interI <- interlist[[i]]
            if(!is.interact(interI))
              stop("Internal error: interlist entry is not an interaction")
            projI <- interI$project
            if(is.null(projI))
              stop(paste("Projection is not yet implemented for a",
                         interI$name))
            p <- projI(coefsI, interI)
            #' p can be NULL (indicating no projection required for interI)
            #' or a single interaction or a list of interactions.
            if(is.null(p)) {
              if(Ncif == 1) return(NULL) # no projection required
              p <- list(NULL)
              nproj[i] <- 0
            } else if(is.interact(p)) {
              p <- list(p)
              nproj[i] <- 1
            } else if(is.list(p) && all(unlist(lapply(p, is.interact)))) {
              nproj[i] <- length(p)
            } else
              stop("Internal error: result of projection had wrong format")
            projlist[[i]] <- p
          }
          #' for interaction i there are nproj[i] **new** interactions to try.
          if(all(nproj == 0))
            return(NULL)
          if(spatstat.options("project.fast")) {
            #' Single interaction required.
            #' Extract first entry from each list
            #' (there should be only one entry, but...)
            qlist <- lapply(projlist, "[[", i=1)
            #' replace NULL entries by corresponding original interactions
            isnul <- unlist(lapply(qlist, is.null))
            if(all(isnul))
              return(NULL)
            if(any(isnul))
              qlist[isnul] <- interlist[isnul]
            names(qlist) <- names(interlist)
            #' build hybrid and return
            result <- do.call("Hybrid", qlist)
            return(result)
          } 
          #' Full case
          result <- list()
          for(i in which(nproj > 0)) {
            ntry <- nproj[i]
            tries <- projlist[[i]]
            for(j in 1:ntry) {
              #' assemble list of component interactions for hybrid
              qlist <- interlist
              qlist[[i]] <- tries[[j]]
              #' eliminate Poisson
              ispois <- unlist(lapply(qlist, is.poisson))
              if(all(ispois)) {
                #' collapse to single Poisson
                h <- Poisson()
              } else {
                if(any(ispois)) qlist <- qlist[!ispois]
                h <- do.call("Hybrid", qlist)
              }
              result <- append(result, list(h))
            }
          }
          #' 'result' is a list of interactions, each a hybrid
          if(length(result) == 1)
            result <- result[[1]]
          return(result)
        },
        irange = function(self, coeffs=NA, epsilon=0, ...) {
          interlist <- self$par
          answer <- 0
          for(i in 1:length(interlist)) {
            interI <- interlist[[i]]
            nameI  <- names(interlist)[[i]]
            nameI. <- paste(nameI, ".", sep="")
            #' find coefficients with prefix that exactly matches nameI.
            if(all(is.na(coeffs)))
              Crelevant <- NA
            else {
              Cname  <- names(coeffs)
              prefixlength <- nchar(nameI.)
              Cprefix <- substr(Cname, 1, prefixlength)
              relevant <- (Cprefix == nameI.)
              #' extract them
              Crelevant <- coeffs[relevant]
              names(Crelevant) <-
                substr(Cname[relevant], prefixlength+1, max(nchar(Cname)))
            }
            #' compute reach 
            reachI <- interI$irange
            if(is.function(reachI)) {
              resultI <- reachI(interI,
                                coeffs=Crelevant, epsilon=epsilon, ...)
              answer <- max(answer, resultI)
            }
          }
          return(answer)
        },
        version=versionstring.spatstat()
        )
    class(out) <- "interact"
    return(out)
  }

  invokeSelfStart <- function(inte, Y) {
    ss <- inte$selfstart
    if(!is.function(ss)) return(inte)
    return(ss(Y, inte))
  }

  Hybrid
})


is.hybrid <- function(x) { UseMethod("is.hybrid") }

is.hybrid.interact <- function(x) {
  return(is.interact(x) && (x$name == "Hybrid interaction"))
}

is.hybrid.ppm <- function(x) {
  return(is.hybrid(as.interact(x)))
}

splitHybridInteraction <- function(coeffs, inte) {
  # For hybrids, $par is a list of the component interactions,
  # but coeffs is a numeric vector. 
  # Split the coefficient vector into the relevant coeffs for each interaction
  interlist <- inte$par
  N <- length(interlist)
  coeflist <- vector(mode="list", length=N)
  for(i in 1:N) {
    interI <- interlist[[i]]
    # forbid hybrids-of-hybrids - these should not occur anyway
    if(interI$name == "Hybrid interaction")
      stop("A hybrid-of-hybrid interactions is not implemented")
    # nameI is the tag that identifies I-th component in hybrid
    nameI  <- names(interlist)[[i]]
    nameI. <- paste(nameI, ".", sep="")
    # find coefficients with prefix that exactly matches nameI.
    Cname  <- names(coeffs)
    prefixlength <- nchar(nameI.)
    Cprefix <- substr(Cname, 1, prefixlength)
    relevant <- (Cprefix == nameI.)
    # extract coefficients
    #   (there may be none, if this interaction is Poisson or an 'offset')
    coeffsI <- coeffs[relevant]
    # remove the prefix so the coefficients are recognisable to interaction
    if(any(relevant)) 
      names(coeffsI) <-
        substr(Cname[relevant], prefixlength+1, max(nchar(Cname)))
    # store
    coeflist[[i]] <- coeffsI
  }
  names(coeflist) <- names(interlist)
  return(list(coeflist=coeflist, interlist=interlist))
}

Hybrid <- intermaker(Hybrid, list(creator="Hybrid",
                                  name="general hybrid Gibbs process",
                                  par=list("..."=42),
                                  parnames=list("any list of interactions")))
