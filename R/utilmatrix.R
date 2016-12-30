#'
#'   utilmatrix.R
#'
#'  Utilities for matrices and arrays
#'
#'       $Revision: 1.1 $ $Date: 2016/12/30 03:22:46 $
#'

matrowsum <- function(x) {
  # was:    x %*% rep.int(1, ncol(x))
  rowSums(x)
}

matcolsum <- function(x) {
  # was:   rep.int(1, nrow(x)) %*% x
  colSums(x)
}
  
matrowany <- function(x) {
  # currently faster than apply(x, 1, any) for logical arrays
  (rowSums(x) > 0)
}

matrowall <- function(x) {
  # currently faster than apply(x, 1, all) for logical arrays
  (rowSums(x) == ncol(x))
}

matcolany <- function(x) {
  # currently faster than apply(x, 2, any) for logical arrays
  (colSums(x) > 0)
}

matcolall <- function(x) {
  # currently faster than apply(x, 2, all) for logical arrays
  (colSums(x) == nrow(x))
}

########
    # hm, this is SLOWER

apply23sum <- function(x) {
  dimx <- dim(x)
  if(length(dimx) != 3)
    stop("x is not a 3D array")
  result <- array(0, dimx[-1])

  nz <- dimx[3]
  for(k in 1:nz) {
    result[,k] <- matcolsum(x[,,k])
  }
  result
}

######################
#
#   matrixsample         subsample or supersample a matrix
#

matrixsample <- function(mat, newdim, phase=c(0,0), scale, na.value=NA) {
  # 'phase+1' is the position of the [1,1] corner of the new matrix
  #  expressed in the coordinates of the old matrix.
  # 'scale' is the size of one step in the new matrix,
  #  expressed in the coordinates of the old matrix.
  # Both 'phase' and 'scale' can take any real value.
  olddim <- dim(mat)
  if(missing(scale)) scale <- (olddim - 1)/(newdim - 1)
  scale <- ensure2vector(scale)
  newdim  <- ensure2vector(newdim)
  newmat <- matrix(na.value, newdim[1], newdim[2])
  newrow <- 1:newdim[1]
  newcol <- 1:newdim[2]
  oldrow <- round(1 + phase[1] + (newrow-1) * scale[1])
  oldcol <- round(1 + phase[2] + (newcol-1) * scale[2])
  oldrow.ok <- (oldrow >= 1) & (oldrow <= olddim[1])
  oldcol.ok <- (oldcol >= 1) & (oldcol <= olddim[2])
  newmat[oldrow.ok, oldcol.ok] <- mat[oldrow[oldrow.ok],
                                      oldcol[oldcol.ok]]
  return(newmat)
}

# wrangle data.frames

findfirstfactor <- function(x) {
  if(!inherits(x, c("data.frame", "hyperframe")))
    stop("x should be a data frame or hyperframe")
  isfac <- unlist(lapply(as.list(x), is.factor))
  if(!any(isfac)) 
    return(NULL)
  min(which(isfac))
}

firstfactor <- function(x) {
  j <- findfirstfactor(x)
  if(is.null(j)) return(NULL)
  return(x[, j, drop=TRUE])
}

assignDFcolumn <- function(x, name, value, ...) {    # for use in mapply 
  dx <- list(value)
  names(dx) <- name
  data.frame(append(c(as.list(x), dx), list(...)))
}

blockdiagmatrix <- function(...) {
  x <- list(...)
  if(!all(unlist(lapply(x, is.matrix))))
    stop("Some of the arguments are not matrices", call.=FALSE)
  nr <- unlist(lapply(x, nrow))
  nc <- unlist(lapply(x, ncol))
  result <- matrix(0, sum(nr), sum(nc))
  rownames(result) <- unlist(lapply(x, rownames))
  colnames(result) <- unlist(lapply(x, colnames))
  rowend <- cumsum(nr)
  rowstart <- c(0, rowend) + 1
  colend <- cumsum(nc)
  colstart <- c(0, colend) + 1
  for(i in seq_along(x))
    result[ (rowstart[i]):(rowend[i]) , (colstart[i]):(colend[i])] <- x[[i]]
  return(result)
}

blockdiagarray <- function(...) {
  x <- list(...)
  if(!all(unlist(lapply(x, is.array))))
    stop("Some of the arguments are not arrays", call.=FALSE)
  dims <- lapply(x, dim)
  dims1 <- unlist(lapply(dims, "[", i=1))
  if(length(dim1 <- unique(dims1)) > 1)
    stop("Arrays have different extents in first dimension")
  dims2 <- unlist(lapply(dims, "[", i=2))
  dims3 <- unlist(lapply(dims, "[", i=3))
  result <- array(0, dim=c(dim1, sum(dims2), sum(dims3)))
  dn <- lapply(x, dimnames)
  dimnames(result)[[2]] <- unlist(lapply(dn, "[[", i=2))
  dimnames(result)[[3]] <- unlist(lapply(dn, "[[", i=3))
  rowend <- cumsum(dims2)
  rowstart <- c(0, rowend) + 1
  colend <- cumsum(dims3)
  colstart <- c(0, colend) + 1
  for(i in seq_along(x))
    result[ , (rowstart[i]):(rowend[i]) , (colstart[i]):(colend[i])] <- x[[i]]
  return(result)
}

asNumericMatrix <- function(x) {
  ## workaround for strange artefact of as.matrix.data.frame
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x
}


exceedsMaxArraySize <- function(...) {
  (prod(as.numeric(c(...))) > .Machine$integer.max)
}



indexCartesian <- function(nn) {
  # enumerate the elements of the Cartesian product of sets,
  # where nn[i] is the size of the i-th set
  as.matrix(do.call(expand.grid, lapply(nn, seq_len)))
}



ensure3Darray <- function(x) {
  nd <- length(dim(x))
  if(nd == 0) {
    x <- array(x, dim=c(length(x), 1, 1))
  } else if(nd == 2) {
    x <- array(x, dim=c(dim(x), 1))
  } else if(nd > 3) {
    laterdims <- dim(x)[-(1:3)]
    if(any(laterdims != 1))
      stop("Higher-dimensional array cannot be reduced to 3 dimensions")
    x <- array(x, dim=dim(x)[1:3])
  }
  return(x)
}

