#'
#'  treebranches.R
#'
#'  Label branches in a tree
#'
#'  $Revision: 1.2 $ $Date: 2015/05/20 09:36:16 $

#' compute branch labels for each *vertex* in the tree L

treebranchlabels <- local({

  treebranchlabels <- function(L, root=1) {
    stopifnot(inherits(L, "linnet"))
    stopifnot(length(root) == 1)
    V <- L$vertices
    #'    M <- L$m
    #' assign label to each vertex
    e <- rep(NA_character_, npoints(V))
    #' do root
    e[root] <- ""
    #' recurse
    descendtree(L, root, e)
  }

  descendtree <- function(L, at, labels, verbose=FALSE) {
    if(verbose)
      cat(paste("Descending from node", at, "\n"))
    below <- which(L$m[at, ] & is.na(labels))
    while(length(below) == 1) {
      if(verbose)
        cat(paste("Line from", at, paren(labels[at]),
                  "to", below, "\n"))
      labels[below] <- labels[at]
      at <- below
      below <- which(L$m[at, ] & is.na(labels))
    }
    if(length(below) == 0) {
      if(verbose) cat("*\n")
      return(labels)
    }
    if(verbose)
      cat(paste("split into", length(below), "\n"))
    if(length(below) > 26)
      stop("Oops - degree > 27")
    labels[below] <- paste(labels[at], letters[1:length(below)], sep="")
    for(b in below)
      labels <- descendtree(L, b, labels)
    return(labels)
  }

  treebranchlabels
})


#' Function which will return the branch label associated with
#' any point on the network

branchlabelfun <- function(L, root=1) {
  L <- as.linnet(L)
  vertexLabels <- treebranchlabels(L, root=root)
  labfrom <- vertexLabels[L$from]
  labto   <- vertexLabels[L$to]
  segmentLabels <- ifelse(nchar(labfrom) < nchar(labto), labto, labfrom)
  f <- function(x, y, seg, tp) { segmentLabels[seg] }
  fL <- linfun(f, L)
  return(fL)
}

#' convenience function for use in model formulae

begins <- function(x, firstbit) {
  stopifnot(is.character(firstbit) && length(firstbit) == 1)
  n <- nchar(firstbit)
  if(n == 0) rep(TRUE, length(x)) else (substr(x, 1, n) == firstbit)
}

#' extract the sub-tree for a particular label
#' e.g. extractbranch(L, "a") extracts everything whose label begins with 'a'

extractbranch <- function(X, ...) {
  UseMethod("extractbranch")
}

extractbranch.linnet <- function(X, code, labels, ..., which=NULL) {
  L <- X
  V <- L$vertices
  if(!is.null(which)) {
    stopifnot(is.logical(which))
    if(length(which) != npoints(V))
      stop("Argument 'which' is the wrong length")
    vin <- which
  } else {
    if(length(labels) != npoints(V))
      stop("labels vector is the wrong length")
    #' which vertices are included
    #'    (a) vertices with the right initial code
    vin <- (substr(labels, 1, nchar(code)) == code)
    #'    (b) the apex
    isneighbour <- (rowSums(L$m[,vin]) > 0)
    apexcode <- if(nchar(code) > 1) substr(code, 1, nchar(code)-1) else ""
    vin <- vin | (isneighbour & (labels == apexcode))
  }
  #' which edges are included
  ein <- vin[L$from] & vin[L$to]
  #' new serial numbers for vertices
  vId <- cumsum(vin)
  #' pack up
  sparse <- L$sparse
  out <- list(vertices=V[vin],
              m=L$m[vin,vin],
              lines=L$lines[ein],
              from=vId[L$from[ein]], to=vId[L$to[ein]],
              dpath=if(sparse) NULL else L$dpath[vin,vin],
              sparse=sparse,
              window=V$window)
  class(out) <- c("linnet", class(out))
  #' pre-compute circumradius
  if(sparse)
    out$circumradius <- circumradius(out)
  attr(out, "which") <- vin
  return(out)  
}

extractbranch.lpp <- function(X, code, labels, ..., which=NULL) {
  L <- as.linnet(X)
  #' make sub-network
  if(missing(code)) code <- NULL
  if(missing(labels)) labels <- NULL
  Lnew <- extractbranch(L, code, labels, which=which)
  #' which vertices are included
  vin <- attr(Lnew, "vin")
  #' which edges are included
  ein <- vin[L$from] & vin[L$to]
  #' which data points are included
  xin <- ein[coords(X)$seg]
  #' new serial numbers for edges
  eId <- cumsum(ein)
  #' construct subset
  Xnew <- X[xin]
  Xnew$domain <- Lnew
  #' apply new serial numbers to segment map
  coords(Xnew)$seg <- eId[coords(Xnew)$seg]
  #'
  return(Xnew)  
}

deletebranch <- function(X, ...) {
  UseMethod("deletebranch")
}

deletebranch.linnet <- function(X, code, labels, ...) {
  L <- X
  V <- L$vertices
  if(length(labels) != npoints(V))
    stop("labels vector is the wrong length")
  #' which vertices are retained
  vkeep <- (substr(labels, 1, nchar(code)) != code)
  #' which edges are retained
  ekeep <- vkeep[L$from] & vkeep[L$to]
  #' new serial numbers for vertices
  vId <- cumsum(vkeep)
  #' pack up
  sparse <- L$sparse
  out <- list(vertices=V[vkeep],
              m=L$m[vkeep,vkeep],
              lines=L$lines[ekeep],
              from=vId[L$from[ekeep]], to=vId[L$to[ekeep]],
              dpath=if(sparse) NULL else L$dpath[vkeep,vkeep],
              sparse=sparse,
              window=V$window)
  class(out) <- c("linnet", class(out))
  #' recompute circumradius
  if(sparse)
    out$circumradius <- circumradius(out)
  attr(out, "which") <- vkeep
  return(out)  
}

deletebranch.lpp <- function(X, code, labels, ...) {
  #' make sub-network
  L <- as.linnet(X)
  Lnew <- deletebranch(L, code=code, labels=labels)
  #' which vertices are retained
  vkeep <- attr(Lnew, "which")
  #' which edges are retained
  ekeep <- vkeep[L$from] & vkeep[L$to]
  #' which data points are retained
  xin <- ekeep[coords(X)$seg]
  #' new serial numbers for vertices
  #        vId <- cumsum(vkeep)
  #' new serial numbers for edges
  eId <- cumsum(ekeep)
  #' construct subset
  Xnew <- X[xin]
  Xnew$domain <- Lnew
  #' apply new serial numbers to segment map
  coords(Xnew)$seg <- eId[coords(Xnew)$seg]
  #'
  return(Xnew)  
}

treeprune <- function(X, root=1, level=0){
  ## collect names of branches to be pruned
  tb <- treebranchlabels(as.linnet(X), root=root)
  keep <- (nchar(tb) <= level)
  Y <- extractbranch(X, which=keep)
  return(Y)
}

