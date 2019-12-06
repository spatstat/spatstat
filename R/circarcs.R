#'
#' circarcs.R
#'
#'  Circular Arcs
#'
#'   An interval on the circle is specified by [left, right]
#'   meaning the arc starting at 'left' going anticlockwise until 'right'.
#'   Here 'left' and 'right' are angles in radians (mod 2*pi)
#'   from the positive x-axis.
#' 
#'   $Revision: 1.5 $ $Date: 2019/12/06 06:15:29 $

check.arc <- function(arc, fatal=TRUE) {
  if(is.numeric(arc) && length(arc) == 2)
    return(TRUE)
  if(fatal)
    stop("arc should be a numeric vector of length 2")
  return(FALSE)
}

inside.arc <- function(theta, arc) {
  check.arc(arc)
  arc <- arc %% (2*pi)
  theta <- theta %% (2*pi)
  if(arc[1] <= arc[2]) {
    #' arc does not cross the positive x-axis
    result <- (arc[1] <= theta) & (theta <= arc[2])
  } else {
    #' arc crosses the positive x-axis
    result <- (arc[1] <= theta) | (theta <= arc[2])
  }
  return(result)
}

circunion <- function(arcs) {
  stopifnot(is.list(arcs))
  nothing <- list()
  everything <- list(c(0, 2*pi))
  if(length(arcs) == 0) return(nothing)
  lapply(arcs, check.arc)
  #' extract all endpoints
  allends <- unlist(arcs)
  allends <- sortunique(as.numeric(allends) %% (2*pi))
  #' compute midpoints between each successive pair of endpoints (mod 2pi)
  midpts <- allends + diff(c(allends, allends[1] + 2*pi))/2
  #' determine which midpoints lie inside one of the arcs
  midinside <- Reduce("|", lapply(arcs, inside.arc, theta=midpts))
  zeroinside <- any(sapply(arcs, inside.arc, theta=0))
  if(!any(midinside) && !zeroinside)
    return(nothing)
  if(all(midinside) && zeroinside)
    return(everything)
  result <- nothing
  if(zeroinside) {
    #' First deal with the connected component containing 0
    #' Scan clockwise from 2*pi for left endpoint of interval
    n <- length(midinside)
    ileft <- (max(which(!midinside)) %% n) + 1L
    aleft <- allends[ileft]
    #' then anticlockwise for right endpoint
    iright <- min(which(!midinside))
    aright <- allends[iright]
    #' save this interval
    result <- append(result, list(c(aleft, aright)))
    #' remove data from consideration
    seqn <- seq_len(n)
    retain <- seqn > iright & seqn < (ileft-1L)
    midinside <- midinside[retain]
    allends   <- allends[retain]
  }
  #' Now scan anticlockwise for first midpoint that is inside the union
  while(any(midinside)) {
    ileft <- min(which(midinside))
    toright  <- (seq_along(midinside) > ileft)
    iright <- min(c(length(allends), which(!midinside & toright)))
    aleft <- allends[ileft]
    aright <- allends[iright]
    #' save this interval
    result <- append(result, list(c(aleft, aright)))
    #' throw away points that are not endpoints of the union
    midinside <- midinside[seq_along(midinside) > iright]
    allends   <- allends[seq_along(allends) > iright]
  }
  return(result)
}

# plotarc <- function(arc, ..., add=TRUE, lwd=3, rad=1){
#   if(!add || is.null(dev.list()))
#     plot(disc(), main="")
#   if(diff(arc) < 0)
#     arc[2] <- arc[2] + 2*pi
#   ang <- seq(arc[1], arc[2], by=0.01)
#   lines(rad * cos(ang), rad * sin(ang), ..., lwd=lwd)
# }
# 
# plotarcs <- function(arcs, ..., rad=1, jitter=FALSE, add=FALSE) {
#   if(length(rad) == 1) rad <- rep(rad, length(arcs))
#   if(jitter) rad <- rad * seq(0.9, 1.05, length=length(rad))
#   rad <- as.list(rad)
#   if(!add) plot(disc(), main="")
#   mapply(plotarc, arc=arcs,rad=rad, MoreArgs=list(...))
#   invisible(NULL)
# }
#   
# runifarc <- function(n=1, maxlen=pi) {
#   replicate(n, runif(1, 0, 2*pi) + c(0, runif(1, 0, maxlen)), simplify=FALSE)
# }
# 
# tryit <- function(n=5, maxlen=pi) {
#   a <- runifarc(n, maxlen=maxlen)
#   plotarcs(circunion(a), col=3, jitter=FALSE, lwd=6)
#   plotarcs(a, jitter=TRUE, lwd=2, add=TRUE)
# }
