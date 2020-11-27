#'
#'  logiquad.R
#'
#'  Quadrature schemes for logistic method
#'
#'  $Revision: 1.1 $  $Date: 2020/11/27 03:02:33 $

logi.dummy <- function(X, dummytype = "stratrand", nd = NULL, mark.repeat = FALSE, ...){
  ## Resolving nd inspired by default.n.tiling
  if(is.null(nd)){
    nd <- spatstat.options("ndummy.min")
    if(inherits(X, "ppp"))
      nd <- pmax(nd, 10 * ceiling(2 * sqrt(X$n)/10))
  }
  nd <- ensure2vector(nd)
  marx <- is.multitype(X)
  if(marx)
    lev <- levels(marks(X))
  if(marx && mark.repeat){
    N <- length(lev)
    Dlist <- inDlist <- vector("list", N)
  } else{
    N <- 1
  }
  W <- as.owin(X)
  type <- match.arg(dummytype, c("stratrand", "binomial", "poisson", "grid", "transgrid"))
  B <- boundingbox(W)
  rho <- nd[1L]*nd[2L]/area(B)
  Dinfo <- list(nd=nd, rho=rho, how=type)
  ## Repeating dummy process for each mark type 1:N (only once if unmarked or mark.repeat = FALSE)
  for(i in 1:N){
    switch(type,
           stratrand={
             D <- as.ppp(stratrand(B, nd[1L], nd[2L]), W = B)
             inD <- which(inside.owin(D, w = W))
             D <- D[W]
             inD <- paste(i,inD,sep="_")
           },
           binomial={
             D <- runifpoint(nd[1L]*nd[2L], win=B)
             D <- D[W]
           },
           poisson={
             D <- rpoispp(rho, win = W)
           },
           grid={
             D <- as.ppp(gridcenters(B, nd[1L], nd[2L]), W = B)
             inD <- which(inside.owin(D, w = W))
             D <- D[W]
             inD <- paste(i,inD,sep="_")
           },
           transgrid={
             D <- as.ppp(gridcenters(B, nd[1L], nd[2L]), W = B)
             dxy <- c(diff(D$window$xrange),diff(D$window$yrange))/(2*nd)
             coords(D) <- coords(D)+matrix(runif(2,-dxy,dxy),npoints(D),2,byrow=TRUE)
             inD <- which(inside.owin(D, w = W))
             D <- D[W]
             inD <- paste(i,inD,sep="_")
           },
         stop("unknown dummy type"))
    if(marx && mark.repeat){
      marks(D) <- factor(lev[i], levels = lev)
      Dlist[[i]] <- D
      if(type %in% c("stratrand","grid","transgrid"))
        inDlist[[i]] <- inD
    }
  }
  if(marx && mark.repeat){
    inD <- Reduce(append, inDlist)
    D <- Reduce(superimpose, Dlist)
  }
  if(type %in% c("stratrand","grid","transgrid"))
    Dinfo <- append(Dinfo, list(inD=inD))
  if(marx && !mark.repeat){
    marks(D) <- sample(factor(lev, levels=lev), npoints(D), replace = TRUE)
    Dinfo$rho <- Dinfo$rho/length(lev)
  }
  attr(D, "dummy.parameters") <- Dinfo
  return(D)
}

quadscheme.logi <- function(data, dummy, dummytype = "stratrand", nd = NULL, mark.repeat = FALSE, ...){
  data <- as.ppp(data)
  ## If dummy is missing we generate dummy pattern with logi.dummy.
  if(missing(dummy))
    dummy <- logi.dummy(data, dummytype, nd, mark.repeat, ...)
  Dinfo <- attr(dummy, "dummy.parameters")
  D <- as.ppp(dummy)
  if(is.null(Dinfo))
    Dinfo <- list(how="given", rho=npoints(D)/(area(D)*markspace.integral(D)))
  ## Weights:
  n <- npoints(data)+npoints(D)
  w <- area(Window(data))/n
  Q <- quad(data, D, rep(w,n), param=Dinfo)
  class(Q) <- c("logiquad", class(Q))
  return(Q)
}

summary.logiquad <- function(object, ..., checkdup=FALSE) {
  verifyclass(object, "logiquad")
  s <- list(
       data  = summary.ppp(object$data, checkdup=checkdup),
       dummy = summary.ppp(object$dummy, checkdup=checkdup),
       param = object$param)
  class(s) <- "summary.logiquad"
  return(s)
}

print.summary.logiquad <- function(x, ..., dp=3) {
  cat("Quadrature scheme (logistic) = data + dummy\n")
  Dinfo <- x$param
  if(is.null(Dinfo))
    cat("created by an unknown function.\n")
  cat("Data pattern:\n")
  print(x$data, dp=dp)

  cat("\n\nDummy pattern:\n")
  # How they were computed
    switch(Dinfo$how,
           stratrand={
             cat(paste("(Stratified random dummy points,",
                       paste(Dinfo$nd, collapse=" x "),
                       "grid of cells)\n"))
           },
           binomial={
             cat("(Binomial dummy points)\n")
           },
           poisson={
             cat("(Poisson dummy points)\n")
           },
           grid={
             cat(paste("(Fixed grid of dummy points,",
                       paste(Dinfo$nd, collapse=" x "),
                       "grid)\n"))
           },
           transgrid={
             cat(paste("(Random translation of fixed grid of dummy points,",
                       paste(Dinfo$nd, collapse=" x "),
                       "grid)\n"))
           },
           given=cat("(Dummy points given by user)\n")
       )
  # Description of them
  print(x$dummy, dp=dp)

  return(invisible(NULL))
}

