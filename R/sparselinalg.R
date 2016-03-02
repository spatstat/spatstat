#'
#'    sparselinalg.R
#'
#'   Counterpart of linalg.R for sparse matrices/arrays
#'
#'   Functions work for 'sparseSlab' or 'sparse3Darray'
#' 
#'   $Revision: 1.3 $  $Date: 2016/03/01 09:13:09 $


tensor1x1 <- function(A, B) {
  ## equivalent of tensor(A, B, 1, 1)
  ## when A is a vector and B is a sparse array.
  stopifnot(length(dim(B)) == 3)
  A <- as.vector(as.matrix(A))
  stopifnot(length(A) == dim(B)[1])
  if(is.array(B)) {
    result <- tensor(A,B,1,1)
  } else if(inherits(B, "sparse3Darray")) {
    result <- sparseMatrix(i=B$j,
                           j=B$k,
                           x=B$x * A[B$i], # values for same (i,j) are summed
                           dims=dim(B)[-1],
                           dimnames=dimnames(B)[2:3])
  } else if(inherits(B, "sparseSlab")) {
    stopifnot(B$stackdim == 1)
    result <- Reduce("+", mapply("*", as.list(A), B$m))
  } else stop("Format of B not understood", call.=FALSE)
  return(result)
}

sumsymouterSparse <- function(x, w=NULL, dbg=FALSE) {
  dimx <- dim(x)
  if(length(dimx) != 3) stop("x should be a 3D array")
  stopifnot(dim(x)[2] == dim(x)[3])
  if(!is.null(w)) {
    stopifnot(inherits(w, "sparseMatrix"))
    stopifnot(all(dim(w) == dim(x)[2:3]))
  }
  m <- dimx[1]
  n <- dimx[2]
  if(inherits(x, "sparse3Darray")) {
    df <- data.frame(i = x$i - 1L,  # need 0-based indices
                     j = x$j - 1L,
                     k = x$k - 1L,
                     value = x$x)
  } else if(inherits(x, "sparseSlab")) {
    if(x$stackdim != 1)
      stop("x should be stacked along the first dimension", call.=FALSE)
    #' Extract triplet representation of x
    mlist <- lapply(x$m, as, Class="TsparseMatrix")
    dflist <- list()
    for(i in seq_along(mlist)) {
      mi <- mlist[[i]]
      dflist[[i]] <- data.frame(i=rep(i-1, length(mi@i)),
                                j=mi@i, k=mi@j, value=mi@x)
      #' note these indices are 0-based
    }
    df <- Reduce(rbind, dflist)
  } else stop("x is not a recognised kind of sparse array")
  # trivial?
  if(nrow(df) < 2) {
    y <- matrix(0, m, m)
    dimnames(y) <- rep(dimnames(x)[1], 2)
    return(y)
  }
  # order by increasing j, then k
  oo <- with(df, order(j, k, i))
  df <- df[oo, ]
  # now provide ordering by increasing k then j
  ff <- with(df, order(k,j,i))
  #
  if(dbg) {
    cat("----------------- Data ---------------------\n")
    print(df)
    cat("-------------- Reordered data --------------\n")
    print(df[ff,])
    cat("Calling......\n")
  }
  if(is.null(w)) {
    z <- .C("CspaSumSymOut",
            m = as.integer(m),
            n = as.integer(n),
            lenx = as.integer(nrow(df)),
            ix = as.integer(df$i), # indices are already 0-based
            jx = as.integer(df$j),
            kx = as.integer(df$k),
            x  = as.double(df$value),
            flip = as.integer(ff - 1L), # convert 1-based to 0-based
            y  = as.double(numeric(m * m)))
  } else {
    # extract triplet representation of w
    w <- as(w, Class="TsparseMatrix")
    dfw <- data.frame(j=w@i, k=w@j, w=w@x)
    woo <- with(dfw, order(j, k))
    dfw <- dfw[woo, , drop=FALSE]
    z <- .C("CspaWtSumSymOut",
            m = as.integer(m),
            n = as.integer(n),
            lenx = as.integer(nrow(df)),
            ix = as.integer(df$i), # indices are already 0-based
            jx = as.integer(df$j),
            kx = as.integer(df$k),
            x  = as.double(df$value),
            flip = as.integer(ff - 1L),  # convert 1-based to 0-based
            lenw = as.integer(nrow(dfw)),
            jw = as.integer(dfw$j),
            kw = as.integer(dfw$k),
            w = as.double(dfw$w),
            y  = as.double(numeric(m * m)))
  }
  y <- matrix(z$y, m, m)
  dimnames(y) <- rep(dimnames(x)[1], 2)
  return(y)
}

