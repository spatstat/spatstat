#
# replace.ppp.R
#


"[<-.ppp" <-
  function(x, i, j, value) {
    verifyclass(x, "ppp")
    verifyclass(value, "ppp")
    
    if(missing(i) && missing(j))
      return(value)

    if(missing(i)) {
      message("The use of argument j in [<-.ppp is deprecated; use argument i")
      # invoke code below
      x[j] <- value
      return(x)
    }

    xmf <- markformat(x)
    vmf <- markformat(value)
    if(xmf != vmf) {
      if(xmf == "none")
        stop("Replacement points are marked, but x is not marked")
      else if(vmf == "none")
        stop("Replacement points have no marks, but x is marked")
      else
        stop("Format of marks in replacement is incompatible with original")
    }
    
    if(inherits(i, "owin")) {
      win <- i
      vok <- inside.owin(value$x, value$y, win)
      if(!all(vok)) {
        warning("Replacement points outside the specified window were deleted")
        value <- value[vok]
      }
      # convert to vector index
      i <- inside.owin(x$x, x$y, win)
    }
    if(!is.vector(i))
      stop("Unrecognised format for subset index i")
    
    # vector index
    # determine index subset
    n <- x$n
    SUB <- seq_len(n)[i]
    # anything to replace?
    if(length(SUB) == 0)
      return(x)
    # sanity checks
    if(any(is.na(SUB)))
      stop("Invalid subset: the resulting subscripts include NAs")
    # exact replacement of this subset?
    if(value$n == length(SUB)) {
      x$x[SUB] <- value$x
      x$y[SUB] <- value$y
      switch(xmf,
             none={},
             list=,
             vector={ x$marks[SUB] <- value$marks },
             dataframe={ x$marks[SUB,] <- value$marks })
    } else 
      x <- superimpose(x[-SUB], value, W=x$window)

    if(!missing(j)) {
      warning("The use of argument j in [<-.ppp is deprecated; use argument i")
      # invoke code above
      x[j] <- value
    }
      
    return(x)
}
