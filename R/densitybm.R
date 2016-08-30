#'
#'   densitybm.R
#'
#'   Barry-McIntyre density estimator
#'
#'   $Revision: 1.2 $ $Date: 2016/08/30 09:38:02 $

densitybm <- function(X, sigma, eps=NULL) {
  stopifnot(is.ppp(X))
  check.1.real(sigma)
  if(is.null(eps)) { eps <- shortside(Frame(X))/128 } else check.1.real(eps)
  Y <- pixellate(X, eps=eps)
  if(with(Y, abs(ystep/xstep - 1) > 0.01))
    stop("Internal error: unequal x and y step lengths")
  Y <- Y/with(Y, xstep * ystep)
  v <- as.matrix(Y)
  #' construct adjacency matrix
  m <- gridadjacencymatrix(dim(v))
  u <- as.vector(v)
  #' restrict to window
  if(anyNA(u)) {
    ok <- !is.na(u)
    m <- m[ok,ok,drop=FALSE]
    u <- u[ok]
  }
  #' construct iteration matrix
  degree <- colSums(m)
  alpha <- 1/5
  A <- alpha * m
  diag(A) <- 1 - alpha * degree
  #' compute variance of one iteration
  stepvar <- 4 * alpha * with(Y, mean(c(xstep,ystep)^2))
  #' number of steps required
  nsteps <- ceiling(sigma^2/stepvar)
  #' run
  U <- u
  for(istep in 1:nsteps) U <- A %*% U
  #' pack up
  Z <- Y
  Z[] <- as.vector(U)
  return(Z)
}


gridadjacencymatrix <- function(dims) {
  dims <- ensure2vector(dims)
  nr <- dims[1]
  nc <- dims[2]
  n <- prod(dims)
  serial <- matrix(1:n, nr, nc)
  allbutlastrow  <- as.vector(serial[-nr,    , drop=FALSE])
  allbutfirstrow <- as.vector(serial[ -1,    , drop=FALSE])
  allbutlastcol  <- as.vector(serial[   , -nc, drop=FALSE])
  allbutfirstcol <- as.vector(serial[   ,  -1, drop=FALSE])
  m <- sparseMatrix(i=integer(0), j=integer(0), x=logical(0), dims=c(n,n))
  m[cbind(allbutfirstrow, allbutlastrow)] <- TRUE
  m[cbind(allbutfirstcol, allbutlastcol)] <- TRUE
  m <- m | t(m)
  return(m)
}


