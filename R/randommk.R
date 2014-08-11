#
#
#   randommk.R
#
#   Random generators for MULTITYPE point processes
#
#   $Revision: 1.30 $   $Date: 2012/08/20 03:48:56 $
#
#   rmpoispp()   random marked Poisson pp
#   rmpoint()    n independent random marked points
#   rmpoint.I.allim()  ... internal
#   rpoint.multi()   temporary wrapper 
#
"rmpoispp" <-
  function(lambda, lmax=NULL, win = owin(c(0,1),c(0,1)),
           types, ...) {
    # arguments:
    #     lambda  intensity:
    #                constant, function(x,y,m,...), image,
    #                vector, list of function(x,y,...) or list of images
    #
    #     lmax     maximum possible value of lambda
    #                constant, vector, or list
    #
    #     win     default observation window (of class 'owin')
    #
    #     types    possible types for multitype pattern
    #    
    #     ...     extra arguments passed to lambda()
    #

    # Validate arguments
    is.numvector <- function(x) {is.numeric(x) && is.vector(x)}
    is.constant <- function(x) {is.numvector(x) && length(x) == 1}
    checkone <- function(x) {
      if(is.constant(x)) {
        if(x >= 0) return(TRUE) else stop("Intensity is negative!")
      }
      return(is.function(x) || is.im(x))
    }
    single.arg <- checkone(lambda)
    vector.arg <- !single.arg && is.numvector(lambda) 
    list.arg <- !single.arg && is.list(lambda)
    if(! (single.arg || vector.arg || list.arg))
      stop(paste("argument", sQuote("lambda"), "not understood"))
    
    if(list.arg && !all(unlist(lapply(lambda, checkone))))
      stop(paste("Each entry in the list",
                 sQuote("lambda"),
                 "must be either a constant, a function or an image"))
    if(vector.arg && any(lambda < 0))
      stop(paste("Some entries in the vector",
                 sQuote("lambda"), "are negative"))


    # Determine & validate the set of possible types
    if(missing(types)) {
      if(single.arg)
        stop(paste(sQuote("types"), "must be given explicitly when",
                   sQuote("lambda"), "is a constant, a function or an image"))
      else
        types <- seq_along(lambda)
    } 

    ntypes <- length(types)
    if(!single.arg && (length(lambda) != ntypes))
      stop(paste("The lengths of", sQuote("lambda"),
                 "and", sQuote("types"), "do not match"))

    factortype <- factor(types, levels=types)

    # Validate `lmax'
    if(! (is.null(lmax) || is.numvector(lmax) || is.list(lmax) ))
      stop(paste(sQuote("lmax"),
                 "should be a constant, a vector, a list or NULL"))
       
    # coerce lmax to a vector, to save confusion
    if(is.null(lmax))
      maxes <- rep(NULL, ntypes)
    else if(is.numvector(lmax) && length(lmax) == 1)
      maxes <- rep(lmax, ntypes)
    else if(length(lmax) != ntypes)
      stop(paste("The length of",
                 sQuote("lmax"),
                 "does not match the number of possible types"))
    else if(is.list(lmax))
      maxes <- unlist(lmax)
    else maxes <- lmax

    # coerce lambda to a list, to save confusion
    lam <- if(single.arg) lapply(1:ntypes, function(x, y){y}, y=lambda)
           else if(vector.arg) as.list(lambda) else lambda

    # Ensure that m can be passed as a single value to function(x,y,m,...)
    slice.fun <- function(x,y,fun,mvalue, ...) {
      m <- if(length(mvalue) == 1) rep(mvalue, length(x)) else mvalue
      result <- fun(x,y,m, ...)
      return(result)
    }
    
    # Simulate
    for(i in 1:ntypes) {
      if(single.arg && is.function(lambda)) 
        # call f(x,y,m, ...)
        Y <- rpoispp(slice.fun, lmax=maxes[i], win=win,
                     fun=lambda, mvalue=types[i], ...)
      else
        # call f(x,y, ...) or use other formats
        Y <- rpoispp(lam[[i]], lmax=maxes[i], win=win, ...)
      Y <- Y %mark% factortype[i]
      X <- if(i == 1) Y else superimpose(X, Y, W=X$window, check=FALSE)
    }

    # Randomly permute, just in case the order is important
    permu <- sample(X$n)
    return(X[permu])
}

# ------------------------------------------------------------------------

"rmpoint" <- function(n, f=1, fmax=NULL, 
                      win = unit.square(), 
                      types, ptypes, ...,
                      giveup = 1000, verbose = FALSE) {
  if(!is.numeric(n))
    stop("n must be a scalar or vector")
  if(any(ceiling(n) != floor(n)))
     stop("n must be an integer or integers")
  if(any(n < 0))
     stop("n must be non-negative")
            
  if(sum(n) == 0) {
    nopoints <- ppp(x=numeric(0), y=numeric(0), window=win, check=FALSE)
    nomarks <- factor(types[numeric(0)], levels=types)
    empty <- nopoints %mark% nomarks
    return(empty)
  }         

  #############
  
  Model <- if(length(n) == 1) {
    if(missing(ptypes)) "I" else "II"
  } else "III"
  
  ##############  Validate f argument
  is.numvector <- function(x) {is.numeric(x) && is.vector(x)}
  is.constant <- function(x) {is.numvector(x) && length(x) == 1}
  checkone <- function(x) {
    if(is.constant(x)) {
      if(x >= 0) return(TRUE) else stop("Intensity is negative!")
    }
    return(is.function(x) || is.im(x))
  }

  single.arg <- checkone(f)
  vector.arg <- !single.arg && is.numvector(f) 
  list.arg <- !single.arg && is.list(f)
  if(! (single.arg || vector.arg || list.arg))
    stop(paste("argument", sQuote("f"), "not understood"))
    
  if(list.arg && !all(unlist(lapply(f, checkone))))
    stop(paste("Each entry in the list", sQuote("f"),
               "must be either a constant, a function or an image"))
  if(vector.arg && any(f < 0))
    stop(paste("Some entries in the vector",
               sQuote("f"), "are negative"))

  # cases where it's known that all types of points 
  # have the same conditional density of location (x,y)
  const.density <- vector.arg ||
             (list.arg && all(unlist(lapply(f, is.constant))))
  same.density <- const.density || (single.arg && !is.function(f))

    
  ################   Determine & validate the set of possible types
  if(missing(types)) {
    if(single.arg && length(n) == 1)
      stop(paste(sQuote("types"), "must be given explicitly when",
                 sQuote("f"), "is a single number, a function or an image and",
                 sQuote("n"), "is a single number"))
    else if(single.arg)
      types <- seq_len(n)
    else 
      types <- seq_along(f)
  }

  ntypes <- length(types)
  if(!single.arg && (length(f) != ntypes))
    stop(paste("The lengths of",
               sQuote("f"), "and", sQuote("types"),
               "do not match"))
  if(length(n) > 1 && ntypes != length(n))
    stop(paste("The lengths of",
               sQuote("n"), "and", sQuote("types"),
               "do not match"))

  factortype <- factor(types, levels=types)
  
  #######################  Validate `fmax'
  if(! (is.null(fmax) || is.numvector(fmax) || is.list(fmax) ))
    stop(paste(sQuote("fmax"),
               "should be a constant, a vector, a list or NULL"))
       
  # coerce fmax to a vector, to save confusion
  if(is.null(fmax))
    maxes <- rep(NULL, ntypes)
  else if(is.constant(fmax))
    maxes <- rep(fmax, ntypes)
  else if(length(fmax) != ntypes)
    stop(paste("The length of", sQuote("fmax"),
               "does not match the number of possible types"))
  else if(is.list(fmax))
    maxes <- unlist(fmax)
  else maxes <- fmax

  # coerce f to a list, to save confusion
  flist <- if(single.arg) lapply(1:ntypes, function(i, f){f}, f=f)
         else if(vector.arg) as.list(f) else f

  #################### START ##################################

  ## special algorithm for Model I when all f[[i]] are images

  if(Model == "I" && !same.density && all(unlist(lapply(flist, is.im))))
    return(rmpoint.I.allim(n, flist, types))

  ## otherwise, first select types, then locations given types
  
  if(Model == "I") {
    # Compute approximate marginal distribution of type
    if(vector.arg)
      ptypes <- f/sum(f)
    else if(list.arg) {
      integratexy <- function(f, win, ...) {
        imag <- as.im(f, win, ...)
        summ <- summary(imag)
        summ$integral
      }
      fintegrals <- unlist(lapply(flist, integratexy, win=win, ...))
      ptypes <- fintegrals/sum(fintegrals)
    } else {
       #single argument
      if(is.constant(f))
        ptypes <- rep(1/ntypes, ntypes)
      else {
        # f is a function (x,y,m)
        # create a counterpart of f that works when m is a single value
        g <- function(xx, yy, ..., m, f) {
          mm <- rep(m, length(xx))
          f(xx, yy, mm, ...)
        }
        # then convert to images and integrate
        fintegrals <-
          unlist(lapply(types,
                        function(typ, ..., win, g) {
                          fim <- as.im(g, W=win, ..., m=typ)
                          summary(fim)$integral
                        },
                        win=win, g=g, f=f))
        # normalise
        ptypes <- fintegrals/sum(fintegrals)
      }
    }
  }

  # Generate marks 

  if(Model == "I" || Model == "II") {
    # i.i.d.: n marks with distribution 'ptypes'
    marques <- sample(factortype, n, prob=ptypes, replace=TRUE)
    nn <- table(marques)
  } else {
    # multinomial: fixed number n[i] of types[i]
    repmarks <- factor(rep(types, n), levels=types)
    marques <- sample(repmarks)
    nn <- n
  }
  ntot <- sum(nn)

  ##############  SIMULATE !!!  #########################

  # If all types have the same conditional density of location,
  # generate the locations using rpoint, and return.
  if(same.density) {
    X <- rpoint(ntot, flist[[1]], maxes[[1]], win=win, ...,
                giveup=giveup, verbose=verbose)
    X <- X %mark% marques
    return(X)
  }
  # Otherwise invoke rpoint() for each type separately
  X <- ppp(numeric(ntot), numeric(ntot), window=win, marks=marques, check=FALSE)

  for(i in 1:ntypes) {
    if(verbose) cat(paste("Type", i, "\n"))
    if(single.arg && is.function(f)) {
      # want to call f(x,y,m, ...)
      # create a counterpart of f that works when m is a single value
      gg <- function(xx, yy, ..., m, fun) {
        mm <- rep(m, length(xx))
        fun(xx, yy, mm, ...)
      }
      Y <- rpoint(nn[i], gg, fmax=maxes[i], win=win,
                  ..., m=factortype[i], fun=f, giveup=giveup, verbose=verbose)
    } else
      # call f(x,y, ...) or use other formats
      Y <- rpoint(nn[i], flist[[i]], fmax=maxes[i], win=win,
                  ..., giveup=giveup, verbose=verbose)
    Y <- Y %mark% factortype[i]
    X[marques == factortype[i]] <- Y
  }
  
  return(X)
}

rmpoint.I.allim <- function(n, f, types) {
  # Internal use only!
  # Generates random marked points (Model I *only*)
  # when all f[[i]] are pixel images.
  #
  # Extract pixel coordinates and probabilities
  get.stuff <- function(imag) {
    w <- as.mask(as.owin(imag))
    dx <- w$xstep
    dy <- w$ystep
    xpix <- as.vector(raster.x(w)[w$m])
    ypix <- as.vector(raster.y(w)[w$m])
    ppix <- as.vector(imag$v[w$m]) # not normalised - OK
    npix <- length(xpix)
    return(list(xpix=xpix, ypix=ypix, ppix=ppix,
                dx=rep(dx,npix), dy=rep(dy, npix),
                npix=npix))
  }
  stuff <- lapply(f, get.stuff)
  # Concatenate into loooong vectors
  xpix <- unlist(lapply(stuff, function(z) { z$xpix }))
  ypix <- unlist(lapply(stuff, function(z) { z$ypix }))
  ppix <- unlist(lapply(stuff, function(z) { z$ppix }))
  dx <- unlist(lapply(stuff, function(z) { z$dx }))
  dy <- unlist(lapply(stuff, function(z) { z$dy }))
  # replicate types
  numpix <- unlist(lapply(stuff, function(z) { z$npix }))
  tpix <- rep(seq_along(types), numpix)
  #
  # sample pixels from union of all images
  #
  npix <- sum(numpix)
  id <- sample(npix, n, replace=TRUE, prob=ppix)
  # get pixel centre coordinates and randomise within pixel
  x <- xpix[id] + (runif(n) - 1/2) * dx[id]
  y <- ypix[id] + (runif(n) - 1/2) * dy[id]
  # compute types
  marx <- factor(types[tpix[id]],levels=types)
  # et voila!
  return(ppp(x, y, window=as.owin(f[[1]]), marks=marx, check=FALSE))
}

#
#     wrapper for Rolf's function
#
rpoint.multi <- function (n, f, fmax=NULL, marks = NULL,
                          win = unit.square(),
                          giveup = 1000, verbose = FALSE,
                          warn=TRUE) {
  no.marks <- is.null(marks) ||
               (is.factor(marks) && length(levels(marks)) == 1)
  if(warn) {
    nhuge <- spatstat.options("huge.npoints")
    if(n > nhuge)
      warning(paste("Attempting to generate", n, "random points"))
  }
  # unmarked case
  if (no.marks) {
    if(is.function(f))
      return(rpoint(n, f, fmax, win, giveup=giveup, verbose=verbose))
    else
      return(rpoint(n, f, fmax, giveup=giveup, verbose=verbose))
  }
  # multitype case
  if(length(marks) != n)
    stop("length of marks vector != n")
  if(!is.factor(marks))
    stop("marks should be a factor")
  types <- levels(marks)
  types <- factor(types, levels=types)
  # generate required number of points of each type
  nums <- table(marks)
  X <- rmpoint(nums, f, fmax, win=win, types=types,
               giveup=giveup, verbose=verbose)
  if(any(table(marks(X)) != nums))
    stop("Internal error: output of rmpoint illegal")
  # reorder them to correspond to the desired 'marks' vector
  Y <- X
  Xmarks <- marks(X)
  for(ty in types) {
    to   <- (marks == ty)
    from <- (Xmarks == ty)
    if(sum(to) != sum(from))
      stop(paste("Internal error: mismatch for mark =", ty))
    if(any(to)) {
      Y$x[to] <- X$x[from]
      Y$y[to] <- X$y[from]
      Y$marks[to] <- ty
    }
  }
  return(Y)
}


  
  

    
