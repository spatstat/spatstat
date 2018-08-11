#
#	marktable.R
#
#	Tabulate mark frequencies in neighbourhood of each point 
#	for multitype point patterns
#
#	$Revision: 1.7 $	$Date: 2015/03/25 03:43:35 $
#
#       Requested by Ian Robertson <igr@stanford.edu>


"marktable" <- 
function(X, R, N, exclude=TRUE, collapse=FALSE) 
{
  verifyclass(X, "ppp")
  if(!is.marked(X, dfok=FALSE))
    stop("point pattern has no marks")
  gotR <- !missing(R) && !is.null(R)
  gotN <- !missing(N) && !is.null(N)
  if(gotN == gotR)
    stop("Exactly one of the arguments N and R should be given")
  stopifnot(is.logical(exclude) && length(exclude) == 1)

  m <- marks(X)
  if(!is.factor(m))
    stop("marks must be a factor")

  if(gotR) {
    stopifnot(is.numeric(R) && length(R) == 1 && R > 0)
    #' identify close pairs
    p <- closepairs(X,R,what="indices")
    pi <- p$i
    pj <- p$j
    if(!exclude) {
      #' add identical pairs
      n <- X$n
      pi <- c(pi, 1:n)
      pj <- c(pj, 1:n)
    }
  } else {
    stopifnot(is.numeric(N) && length(N) == 1)
    ii <- seq_len(npoints(X))
    nn <- nnwhich(X, k=1:N)
    if(N == 1) nn <- matrix(nn, ncol=1)
    if(!exclude)
      nn <- cbind(ii, nn)
    pi <- as.vector(row(nn))
    pj <- as.vector(nn)
  }

  #' tabulate
  if(!collapse) {
    ## table for each point
    i <- factor(pi, levels=seq_len(npoints(X)))
    mj <- m[pj]
    mat <- table(point=i, mark=mj)
  } else {
    #' table by type
    mi <- m[pi]
    mj <- m[pj]
    mat <- table(point=mi, neighbour=mj)
  }
  return(mat)
}

