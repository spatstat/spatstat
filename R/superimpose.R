# superimpose.R
#
# $Revision: 1.28 $ $Date: 2014/07/21 06:40:11 $
#
#
############################# 

superimpose <- function(...) {
  # remove any NULL arguments
  arglist <- list(...)
  if(any(isnull <- sapply(arglist, is.null)))
    return(do.call("superimpose", arglist[!isnull]))
  UseMethod("superimpose")
}

superimpose.default <- function(...) {
  argh <- list(...)
  if(any(sapply(argh, is.lpp)) || any(sapply(argh, inherits, what="linnet")))
    return(superimpose.lpp(...))
  return(superimpose.ppp(...))
}

superimpose.ppp <- function(..., W=NULL, check=TRUE) {
  arglist <- list(...)
  # Check that all "..." arguments have x, y coordinates
  hasxy <- unlist(lapply(arglist, checkfields, L=c("x", "y")))
  if(!all(hasxy)) {
    nbad <- sum(bad <- !hasxy)
    stop(paste(ngettext(nbad, "Argument", "Arguments"),
               commasep(which(bad)),
               ngettext(nbad, "does not", "do not"),
               "have components x and y"))
  }
  
  # concatenate lists of (x,y) coordinates
  XY <- do.call("concatxy", arglist)
  needcheck <- TRUE

  # determine whether there is any window information
  if(!is.owin(W)) {
    # we have to compute the final window
    WXY <- NULL
    Wppp <- NULL
    if(any(isppp <- unlist(lapply(arglist, is.ppp)))) {
      # extract windows from ppp objects
      wins <- unname(lapply(arglist[isppp], as.owin))
      # take union
      Wppp <- if(length(wins) == 1) wins[[1]] else do.call(union.owin, wins)
    } 
    if(is.function(W)) {
      # W is a function like bounding.box.xy or ripras
      # Apply function to the x,y coordinates; it should return an owin
      WXY <- W(XY)
      if(!is.owin(WXY))
        stop("Function W did not return an owin object")
    }
    if(is.character(W)) {
      # character string identifies a function
      pW <- pmatch(W, c("convex", "rectangle", "bbox", "none"))
      if(is.na(pW))
        stop(paste("Unrecognised option W=", sQuote(W)))
      WXY <- switch(pW,
                    convex=ripras(XY),
                    rectangle=ripras(XY, shape="rectangle"),
                    bbox=boundingbox(XY),
                    none=NULL)
      # in these cases we don't need to verify that the points are inside.
      needcheck <- !is.null(WXY)
    }
    if(is.null(WXY) && is.null(Wppp)) {
      # no window information
      return(XY)
    }
    W <- union.owin(WXY, Wppp)
  }
  # extract the marks if any
  nobj <- sapply(arglist, function(x) { length(x$x) })
  marx  <- superimposeMarks(arglist, nobj)
  #
  ppp(XY$x, XY$y, window=W, marks=marx, check=check & needcheck)
}

superimpose.psp <- function(..., W=NULL, check=TRUE) {
  # superimpose any number of line segment patterns
  arglist <- list(...)
  misscheck <- missing(check)

  if(!all(sapply(arglist, is.psp)))
    stop("Patterns to be superimposed must all be psp objects")

  # extract segment coordinates
  matlist <- lapply(lapply(arglist, getElement, name="ends"),
                    asNumericMatrix)
  
  # tack them together
  mat <- do.call("rbind", matlist)

  # determine whether there is any window information
  needcheck <- FALSE
  if(!is.owin(W)) {
    # we have to compute the final window
    WXY <- NULL
#    Wpsp <- NULL
    if(any(ispsp <- unlist(lapply(arglist, is.psp)))) {
      # extract windows from psp objects
      wins <- unname(lapply(arglist[ispsp], as.owin))
      # take union
      Wppp <- if(length(wins) == 1) wins[[1]] else do.call(union.owin, wins)
    }
    if(is.function(W) || is.character(W)) {
      # guess window from x, y coordinates
      XY <- list(x=cbind(mat[,1], mat[,3]),
                 y=cbind(mat[,2], mat[,4]))
      if(is.function(W)) {
        # W is a function like bounding.box.xy or ripras
        # Apply function to the x,y coordinates; it should return an owin
        WXY <- W(XY)
        if(!is.owin(WXY))
          stop("Function W did not return an owin object")
      }
      if(is.character(W)) {
        # character string identifies a function
        pW <- pmatch(W, c("convex", "rectangle", "bbox", "none"))
        if(is.na(pW))
          stop(paste("Unrecognised option W=", sQuote(W)))
        WXY <- switch(pW,
                      convex=ripras(XY),
                      rectangle=ripras(XY, shape="rectangle"),
                      bbox=boundingbox(XY),
                      none=NULL)
      # in these cases we don't need to verify that the points are inside.
        needcheck <- !is.null(WXY)
      }
    }
    W <- union.owin(WXY, Wppp)
  }
  
  # extract marks, if any
  nobj <- sapply(arglist, nsegments)
  marx <- superimposeMarks(arglist, nobj)

  if(misscheck && !needcheck) check <- FALSE
  return(as.psp(mat, window=W, marks=marx, check=check))
}

superimposeMarks <- function(arglist, nobj) {
  # combine marks from the objects in the argument list
  marxlist <- lapply(arglist, marks)
  marx <- do.call(markappend, unname(marxlist))
  nama <- names(arglist)
  if(length(nama) == length(arglist) && all(nzchar(nama))) {
    # arguments are named: use names as (extra) marks
    newmarx <- factor(rep.int(nama, nobj))
    marx <- markcbind(marx, newmarx)
    if(ncol(marx) == 2) {
      ## component marks were not named: call them 'origMarks'
      colnames(marx) <- c("origMarks", "pattern")
    } else colnames(marx)[ncol(marx)] <- "pattern"
  }
  return(marx)
}

#==+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===

  # This function is now deprecated.
superimposePSP <-
  function(..., W=NULL, check=TRUE)
{
  .Deprecated("superimpose","spatstat")
  
  # superimpose any number of line segment patterns
  arglist <- list(...)

  nargue <- length(arglist)
  if(nargue == 0)
    stop("No line segment patterns given")
  
  # catch possible abuses
  if(is.null(W) && any(suspicious <- (names(arglist) == "window"))) {
    id <- min(which(suspicious))
    Win <- arglist[[id]]
    if(is.owin(Win) || is.null(Win)) {
      W <- Win
      arglist <- arglist[-id]
      nargue <- length(arglist)
    }
  }

  # unpack a list
  if(nargue == 1) {
    X <- arglist[[1]]
    if(!inherits(X, "psp") && inherits(X, "list"))
      arglist <- X
  }

  isnull <- unlist(lapply(arglist, is.null))
  arglist <- arglist[!isnull]
  
  if(!all(unlist(lapply(arglist, is.psp))))
    stop("Some of the arguments are not psp objects")
  
  # extract segment coordinates
  matlist <- lapply(arglist, function(x) { as.matrix(x$ends) })
  # tack them together
  mat <- do.call("rbind", matlist)

  # extract marks if any
  marxlist <- lapply(arglist, marks)

  # check on compatibility of marks
  mkfmt <- sapply(marxlist,markformat)
  if(length(unique(mkfmt))>1)
	stop(paste("Marks of some patterns are of different format\n",
                   "  from those of other patterns.\n",sep=""))
  mkfmt <- mkfmt[1]
  if(mkfmt=="dataframe") {
	mcnms <- lapply(marxlist,names)
	cdim  <- sapply(mcnms,length)
	OK    <- length(unique(cdim)) == 1
	if(OK) {
		allInOne <- sapply(mcnms,paste,collapse="")
		OK <- length(unique(allInOne)) == 1
		if(!OK) stop("Data frames of marks have different names.\n")
	} else stop("Data frames of marks have different column dimensions.\n")
  }
 
  # combine the marks
  marx <- switch(mkfmt,
                 none = NULL,
                 vector = {
                   marxlist <- lapply(marxlist,
                                      function(x){as.data.frame.vector(x,nm="v1")})
                   do.call("rbind", marxlist)[,1]
                 },
                 dataframe = do.call("rbind", marxlist))

  # determine window
  if(!is.null(W))
    W <- as.owin(W)
  else {
    # extract windows from psp objects
    Wlist <- lapply(arglist, as.owin)
    # take the union of all the windows
    W <- NULL
    for(i in seq_along(Wlist))
      W <- union.owin(W, Wlist[[i]])
  }

  return(as.psp(mat, window=W, marks=marx, check=check))
}

