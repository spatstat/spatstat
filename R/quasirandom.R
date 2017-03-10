##
##     quasirandom.R
##
##  Quasi-random sequence generators
##
##  $Revision: 1.5 $   $Date: 2014/10/24 00:22:30 $
##

vdCorput <- function(n, base) {
  stopifnot(is.prime(base))
  z <- .C("Corput",
          base=as.integer(base),
          n=as.integer(n),
          result=as.double(numeric(n)),
          PACKAGE = "spatstat")
  return(z$result)
}

Halton <- function(n, bases=c(2,3), raw=FALSE, simplify=TRUE) {
  d <- length(bases)
  if(d==2 && !raw && simplify)
    return(ppp(vdCorput(n, bases[1]),
               vdCorput(n, bases[2]),
               window=owin(), check=FALSE))
  z <- matrix(, nrow=n, ncol=d)
  for(j in 1:d)
    z[,j] <- vdCorput(n, bases[j])
  if(raw || d < 2) return(z)
  b <- do.call(boxx, rep(list(c(0,1)), d))
  return(ppx(z, b, simplify=simplify))
}

Hammersley <- function(n, bases=2, raw=FALSE, simplify=TRUE) {
  d <- length(bases) + 1
  z <- cbind(Halton(n, bases, raw=TRUE), (1:n)/n)
  dimnames(z) <- NULL
  if(raw || d < 2) return(z)
  b <- do.call(boxx, rep(list(c(0,1)), d))
  return(ppx(z, b, simplify=simplify))
}

rQuasi <- function(n, W, type=c("Halton", "Hammersley"), ...) {
  R <- as.rectangle(W)
  type <- match.arg(type)
  X <- switch(type,
              Halton=Halton(n, ...),
              Hammersley=Hammersley(n, ...))
  Y <- ppp(R$xrange[1] + diff(R$xrange) * X$x,
           R$yrange[1] + diff(R$yrange) * X$y,
           window=R, check=FALSE)
  if(!is.rectangle(W))
    Y <- Y[W]
  return(Y)
}

