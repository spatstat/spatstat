#
#
#    concom.R
#
#    $Revision: 1.4 $	$Date: 2016/04/25 02:34:40 $
#
#    The connected component interaction
#
#    Concom()    create an instance of the connected component interaction
#                 [an object of class 'interact']
#	
# -------------------------------------------------------------------
#

Concom <- local({

  connectedlabels <- function(X, R) {
    connected(X, R, internal=TRUE)
  }
  
  countcompo <- function(X, R) {
    length(unique(connectedlabels(X, R)))
  }

   # change in number of components when point i is deleted
  cocoDel <- function(X, R, subset=seq_len(npoints(X))) {
    n <- length(subset)
    ans <- integer(n)
    if(n > 0) {
      cX <- countcompo(X, R)
      for(i in 1:n) 
        ans[i] = countcompo(X[-subset[i]], R) - cX
    }
    return(ans)
  }

  # change in number of components when new point is added

  cocoAdd <- function(U, X, R) {
    U <- as.ppp(U, W=as.owin(X))
    nU <- npoints(U)
    cr <- crosspairs(U, X, R, what="indices")
    lab <- connectedlabels(X, R)
    hitcomp <- tapply(X=lab[cr$j],
                      INDEX=factor(cr$i, levels=1:nU),
                      FUN=unique, 
                      simplify=FALSE)
    nhit <- unname(lengths(hitcomp))
    change <- 1L - nhit
    return(change)
  }

  # connected component potential 
  cocopot <- 
    function(X,U,EqualPairs,pars,correction, ...) {
      bad <- !(correction %in% c("border", "none"))
      if((nbad <- sum(bad)) > 0) 
        warning(paste("The",
                      ngettext(nbad, "correction", "corrections"),
                      commasep(sQuote(correction[!ok])),
                      ngettext(nbad, "is", "are"),
                      "not implemented"))
      n <- U$n
      answer <- numeric(n)
      r <- pars$r
      if(is.null(r)) stop("internal error: r parameter not found")
      dummies <- !(seq_len(n) %in% EqualPairs[,2L])
      if(sum(dummies) > 0)
        answer[dummies] <- -cocoAdd(U[dummies], X, r)
      ii <- EqualPairs[,1L]
      jj <- EqualPairs[,2L]
      answer[jj] <- cocoDel(X, r, subset=ii)
      return(answer + 1L)
    }

  # template object without family, par, version
  BlankCoco <- 
  list(
         name     = "Connected component process",
         creator  = "Concom",
         family   = "inforder.family", # evaluated later
         pot      = cocopot,
         par      = list(r = NULL), # to be filled in
         parnames = "distance threshold",
         init     = function(self) {
                      r <- self$par$r
                      if(!is.numeric(r) || length(r) != 1L || r <= 0)
                       stop("distance threshold r must be a positive number")
                    },
         update = NULL,  # default OK
         print = NULL,    # default OK
         interpret =  function(coeffs, self) {
           logeta <- as.numeric(coeffs[1L])
           eta <- exp(logeta)
           return(list(param=list(eta=eta),
                       inames="interaction parameter eta",
                       printable=signif(eta)))
         },
         valid = function(coeffs, self) {
           eta <- ((self$interpret)(coeffs, self))$param$eta
           return(is.finite(eta))
         },
         project = function(coeffs, self) {
           if((self$valid)(coeffs, self))
             return(NULL)
           return(Poisson())
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           if(anyNA(coeffs))
             return(Inf)
           logeta <- coeffs[1L]
           if(abs(logeta) <= epsilon)
             return(0)
           else
             return(Inf)
         },
       version=NULL # to be added
  )
  class(BlankCoco) <- "interact"

  Concom <- function(r) {
    instantiate.interact(BlankCoco, list(r=r))
  }

  Concom <- intermaker(Concom, BlankCoco)
  
  Concom
})


             
