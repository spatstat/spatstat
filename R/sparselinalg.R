#'
#'    sparselinalg.R
#'
#'   Counterpart of linalg.R for sparse matrices/arrays
#'
#'   $Revision: 1.1 $  $Date: 2016/02/21 09:35:48 $


sumsymouterSparseSlab <- function(x, w=NULL, dbg=FALSE) {
  stopifnot(inherits(x, "sparseSlab"))
  if(x$stackdim != 1)
    stop("x should be stacked along the first dimension", call.=FALSE)
  d <- dim(x)
  m <- d[1]
  n <- d[2]
  if(d[3] != n)
    stop("The second and third dimensions of x should be equal",
         call.=FALSE)
  if(!is.null(w)) {
    stopifnot(inherits(w, "sparseMatrix"))
    stopifnot(all(dim(w) == dim(x)[2:3]))
  }
  # Extract triplet representation of x
  mlist <- lapply(x$m, as, Class="TsparseMatrix")
  dflist <- list()
  for(i in seq_along(mlist)) {
    mi <- mlist[[i]]
    dflist[[i]] <- data.frame(i=rep(i, length(mi@i)),
                              j=mi@i, k=mi@j, value=mi@x)
  }
  df <- Reduce(rbind, dflist)
  # trivial?
  if(nrow(df) < 2)
    return(matrix(0, m, m))
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
            ix = as.integer(df$i - 1L),
            jx = as.integer(df$j - 1L),
            kx = as.integer(df$k - 1L),
            x  = as.double(df$value),
            flip = as.integer(ff - 1L),
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
            ix = as.integer(df$i - 1L),
            jx = as.integer(df$j - 1L),
            kx = as.integer(df$k - 1L),
            x  = as.double(df$value),
            flip = as.integer(ff - 1L),
            lenw = as.integer(nrow(dfw)),
            jw = as.integer(dfw$j - 1L),
            kw = as.integer(dfw$k - 1L),
            w = as.double(dfw$w),
            y  = as.double(numeric(m * m)))
  }
  y <- matrix(z$y, m, m)
  dimnames(y) <- rep(dimnames(x)[1], 2)
  return(y)
}

tensor1x1 <- function(A, B) {
  ## equivalent of tensor(A, B, 1, 1)
  ## when A is a vector and B is a sparseSlab.
  stopifnot(inherits(B, "sparseSlab"))
  stopifnot(length(A) == dim(B)[1])
  stopifnot(B$stackdim == 1)
  z <- Reduce("+", mapply("*", as.list(A), B$m))
  return(z)
}
