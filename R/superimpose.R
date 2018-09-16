# superimpose.R
#
# $Revision: 1.37 $ $Date: 2017/12/30 05:03:15 $
#
#
############################# 

superimpose <- function(...) {
  # remove any NULL arguments
  arglist <- list(...)
  if(any(isnull <- sapply(arglist, is.null)))
    return(do.call(superimpose, arglist[!isnull]))
  UseMethod("superimpose")
}

superimpose.default <- function(...) {
  argh <- list(...)
  #' First expand any arguments which are lists of objects
  argh <- expandSpecialLists(argh, "solist")
  #' Now dispatch
  if(any(sapply(argh, is.lpp)) || any(sapply(argh, inherits, what="linnet")))
    return(do.call(superimpose.lpp, argh))
  if(any(sapply(argh, is.psp)))
    return(do.call(superimpose.psp, argh))
  #' default
  return(do.call(superimpose.ppp, argh))
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
               "have components x and y"),
         call.=FALSE)
  }
  
  # concatenate lists of (x,y) coordinates
  XY <- do.call(concatxy, arglist)
  needcheck <- TRUE

  # determine whether there is any window information
  if(!is.owin(W)) {
    ## we have to compute the final window
    if(is.function(W)) {
      ## W is a function like bounding.box.xy or ripras
      ## Apply function to the x,y coordinates; it should return an owin
      WXY <- W(XY)
      if(!is.owin(WXY))
        stop("Function W did not return an owin object", call.=FALSE)
      W <- WXY
    } else if(is.character(W)) {
      ## character string identifies a function
      pW <- pmatch(W, c("convex", "rectangle", "bbox", "none"))
      if(is.na(pW))
        stop(paste("Unrecognised option W=", sQuote(W)), call.=FALSE)
      WXY <- switch(pW,
                    convex=ripras(XY),
                    rectangle=ripras(XY, shape="rectangle"),
                    bbox=boundingbox(XY),
                    none=NULL)
      # in these cases we don't need to verify that the points are inside.
      needcheck <- !is.null(WXY)
      if(!is.null(WXY)) W <- WXY
    } else if(is.null(W)) {
      if(any(isppp <- unlist(lapply(arglist, is.ppp)))) {
        ## extract windows from ppp objects
        wins <- unname(lapply(arglist[isppp], as.owin))
        ## take union
        W <- if(length(wins) == 1) wins[[1]] else do.call(union.owin, wins)
      } else {
        ## no window information
        return(XY)
      }
    } else stop("Argument W is not understood")
  }
  # extract the marks if any
  nobj <- lengths(lapply(arglist, getElement, name="x"))
  marx  <- superimposeMarks(arglist, nobj)
  #
  ppp(XY$x, XY$y, window=W, marks=marx, check=check & needcheck)
}

superimpose.splitppp <- superimpose.ppplist <-
  function(..., W=NULL, check=TRUE) {
    arglist <- list(...)
    while(any(h <- sapply(arglist, inherits, what=c("splitppp", "ppplist")))) {
      i <- min(which(h))
      arglist <- insertinlist(arglist, i, arglist[[i]])
    }
    do.call(superimpose, append(arglist, list(W=W, check=check)))
  }

superimpose.psp <- function(..., W=NULL, check=TRUE) {
  # superimpose any number of line segment patterns
  arglist <- list(...)
  misscheck <- missing(check)

  if(!all(sapply(arglist, is.psp)))
    stop("Patterns to be superimposed must all be psp objects", call.=FALSE)

  # extract segment coordinates
  matlist <- lapply(lapply(arglist, getElement, name="ends"),
                    asNumericMatrix)
  
  # tack them together
  mat <- do.call(rbind, matlist)

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
          stop("Function W did not return an owin object", call.=FALSE)
      }
      if(is.character(W)) {
        # character string identifies a function
        pW <- pmatch(W, c("convex", "rectangle", "bbox", "none"))
        if(is.na(pW))
          stop(paste("Unrecognised option W=", sQuote(W)), call.=FALSE)
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
    newmarx <- factor(rep.int(nama, nobj), levels=nama)
    marx <- markcbind(marx, newmarx)
    if(ncol(marx) == 2) {
      ## component marks were not named: call them 'origMarks'
      colnames(marx) <- c("origMarks", "pattern")
    } else colnames(marx)[ncol(marx)] <- "pattern"
  }
  return(marx)
}

