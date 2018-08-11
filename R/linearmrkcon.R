#
# linearmrkcon.R
#
# mark connection function & mark equality function for linear networks
#
# $Revision: 1.4 $ $Date: 2017/02/07 08:12:05 $
#

linearmarkconnect <- function(X, i, j, r=NULL, ...) {
  if(!is.multitype(X, dfok=FALSE)) 
	stop("Point pattern must be multitype")
  marx <- marks(X)
  lev <- levels(marx)
  if(missing(i) || is.null(i)) i <- lev[1L] else
    if(!(i %in% lev)) stop(paste("i = ", i , "is not a valid mark"))
  if(missing(j) || is.null(j)) j <- lev[2L] else
    if(!(j %in% lev)) stop(paste("j = ", j , "is not a valid mark"))

  # ensure distance information is present
  X <- as.lpp(X, sparse=FALSE)

  #
  pcfij <- linearpcfcross(X, i, j, r=r, ...)
  pcfall <- linearpcf(X, r=r, ...)

  qi <- mean(marx == i)
  qj <- mean(marx == j)

  result <- eval.fv(qi * qj * pcfij/pcfall)
  
  # rebrand
  result <- rebadge.as.crossfun(result, "p", "L", i, j)
  attr(result, "labl") <- attr(pcfij, "labl")
  return(result)
}

linearmarkequal <- local({
  
  linearmarkequal <- function(X, r=NULL, ...) {
    if(!is.multitype(X, dfok=FALSE)) 
      stop("Point pattern must be multitype")
  
    ## ensure distance information is present
    X <- as.lpp(X, sparse=FALSE)

    lev <- levels(marks(X))
    v <- list()
    for(l in lev) v[[l]] <- linearmarkconnect(X, l, l, r=r, ...)

    result <- Reduce(addfuns, v)
    result <-rebadge.fv(result, 
                        quote(p[L](r)),
                        new.fname=c("p", "L"))
    attr(result, "labl") <- attr(v[[1L]], "labl")
    return(result)
  }

  addfuns <- function(f1, f2) eval.fv(f1 + f2)

  linearmarkequal
})


