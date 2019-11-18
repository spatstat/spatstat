#
# rshift.psp.R
#
#  $Revision: 1.7 $  $Date: 2019/11/18 06:22:50 $
#


rshift.psp <- function(X, ..., group=NULL, which=NULL) {
  verifyclass(X, "psp")
  
  # process arguments
  W <- rescue.rectangle(X$window)
  arglist <- handle.rshift.args(W, ..., edgedefault="erode")
  radius <- arglist$radius
  width  <- arglist$width
  height <- arglist$height
  edge   <- arglist$edge
  clip   <- arglist$clip
  if(W$type != "rectangle")
    stop("Not yet implemented for non-rectangular windows")
  if(edge != "erode")
    stop(paste("Only implemented for edge=", dQuote("erode")))

  # split into groups
  if(is.null(group))
    Y <- list(X)
  else {
    stopifnot(is.factor(group))
    stopifnot(length(group) == X$n)
    Y <- lapply(levels(group),
                function(l, Z, group) {Z[group == l]},
                Z=X, group=group)
  }

  ############ loop ################
  result <- NULL
  
  for(i in seq_along(Y)) {
    
    Z <- Y[[i]]
    
    # generate random translation vector
    if(!is.null(radius)) 
      jump <- runifdisc(1, radius=radius)
    else {
      jump <- list(x=runif(1, min=0, max=width),
                   y=runif(1, min=0, max=height))
    }
    # translate segments
    Zsh <- shift(Z, c(jump$x, jump$y))
    Zsh$window <- W

    # append to result
    result <- append.psp(result, Zsh)
  }

  # clip 
  if(!is.null(clip))
   result <- result[clip]

  return(result)
}

