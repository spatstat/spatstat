#
#   rshift.R
#
#   random shift with optional toroidal boundary
#
#   $Revision: 1.17 $   $Date: 2016/02/11 10:17:12 $
#
#
rshift <- function(X, ...) {
  UseMethod("rshift")
}

rshift.splitppp <- function(X, ..., which=seq_along(X))
{
  verifyclass(X, "splitppp")
  if("group" %in% names(list(...)))
    stop(paste("argument", sQuote("group"),
               "not implemented for splitppp objects"))

  if(is.null(which)) {
    iwhich <- which <- seq_along(X)
  } else {
    id <- seq_along(X)
    names(id) <- names(X)
    iwhich <- id[which]
    if(length(iwhich) == 0)
      stop(paste("Argument", sQuote("which"), "did not match any marks"))
  }
  
  # validate arguments and determine common clipping window
  arglist <- handle.rshift.args(X[[1]]$window, ..., edgedefault="torus")

  if(!is.null(clip <- arglist$clip)) {
    # clip the patterns that are not to be shifted
    if(length(iwhich) < length(X)) 
      X[-iwhich] <- lapply(X[-iwhich], "[.ppp", i=clip)
  }
  # perform shift on selected patterns
  # (setting group = NULL ensures each pattern is not split further)
  shiftXsub <- do.call(lapply, append(list(X[iwhich], rshift.ppp, group=NULL),
                                        arglist))
  # put back
  X[iwhich] <- shiftXsub

  return(X)
}

rshift.ppp <- function(X, ..., which=NULL, group)
{
  verifyclass(X, "ppp")
  
  # validate arguments and determine common clipping window
  arglist <- handle.rshift.args(X$window, ..., edgedefault="torus")

  # default grouping
  #   (NULL is not the default)
  #   (NULL means all points shifted in parallel)
  if(missing(group))
    group <- if(is.multitype(X)) marks(X) else NULL

  # if no grouping, use of `which' is undefined
  if(is.null(group) && !is.null(which))
    stop(paste("Cannot apply argument", sQuote("which"),
               "; no grouping defined"))

  # if grouping, use split
  if(!is.null(group)) {
    Y <- split(X, group)
    split(X, group) <- do.call(rshift.splitppp,
                               append(list(Y, which=which),
                                      arglist))
    return(X)
  } 
    
  # ungrouped point pattern
  # shift all points in parallel

  # recover arguments
  radius <- arglist$radius
  width  <- arglist$width
  height <- arglist$height
  edge   <- arglist$edge
  clip   <- arglist$clip
 
  W <- X$window
  W <- rescue.rectangle(W)
  if(W$type != "rectangle" && edge=="torus")
    stop("Torus (periodic) boundary is only meaningful for rectangular windows")

  # generate random translation vector
  
  if(!is.null(radius)) 
    jump <- runifdisc(1, radius=radius)
  else {
    jump <- list(x=runif(1, min=0, max=width),
                 y=runif(1, min=0, max=height))
  }

  # translate points
  x <- X$x + jump$x
  y <- X$y + jump$y

  # wrap points
  if(edge == "torus") {
    xr <- W$xrange
    yr <- W$yrange
    Wide <- diff(xr)
    High <- diff(yr)
    x <- xr[1] + (x - xr[1]) %% Wide
    y <- yr[1] + (y - yr[1]) %% High
  }

  # put back into point pattern
  X$x <- x
  X$y <- y

  # clip to window
  if(!is.null(clip))
    X <- X[clip]

  return(X)
}


handle.rshift.args <- function(W, ...,
                               radius=NULL, width=NULL, height=NULL,
                               edge=NULL, clip=NULL, edgedefault)
{
  verifyclass(W, "owin")
  W <- rescue.rectangle(W)
  
  if(length(aargh <- list(...)) > 0)
    stop(paste("Unrecognised arguments:",
               paste(names(aargh), collapse=",")))
  
  if(!is.null(radius)) {
    # radial generator
    if(!(is.null(width) && is.null(height)))
    stop(paste(sQuote("radius"), "is incompatible with",
               sQuote("width"), "and", sQuote("height")))
  } else {
    # rectangular generator
    if(is.null(width) != is.null(height))
      stop("Must specify both width and height, if one is specified")
    if(is.null(width)) width <- diff(W$xrange)
    if(is.null(height)) height <- diff(W$yrange)
  }
  
  if(is.null(edge))
    edge <- edgedefault
  else if(!(edge %in% c("torus", "erode", "none")))
    stop(paste("Unrecognised option erode=", sQuote(edge)))

  # determine whether clipping window is needed
  if(is.null(clip))
    clip <- switch(edge,
                   torus= NULL,
                   none= W,
                   erode={
                     if(!is.null(radius))
                       erosion.owin(W, radius)
                     else if(W$type == "rectangle")
                       trim.rectangle(W, width, height)
                     else
                       erosion.owin(W, max(width, height))
                   })

  return(list(radius=radius, width=width, height=height,
              edge=edge, clip=clip))
}

rtoro <- function(X, which=NULL, radius=NULL, width=NULL, height=NULL) {
  .Deprecated("rshift", package="spatstat")
  rshift(X, which=which, radius=radius, width=width, height=height)
}
