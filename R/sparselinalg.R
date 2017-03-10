#'
#'    sparselinalg.R
#'
#'   Counterpart of linalg.R for sparse matrices/arrays
#'
#' 
#'   $Revision: 1.9 $  $Date: 2016/04/25 02:34:40 $

marginSums <- function(X, MARGIN) {
  #' equivalent to apply(X, MARGIN, sum)
  if(length(MARGIN) == 0) return(sum(X))
  if(is.array(X) || is.matrix(X))
    return(apply(X, MARGIN, sum))
  dimX <- dim(X)
  if(length(MARGIN) == length(dimX)) 
    return(aperm(X, MARGIN))
  if(any(huh <- (MARGIN < 0 | MARGIN > length(dimX))))
    stop(paste(commasep(sQuote(paste0("MARGIN=", MARGIN[huh]))),
               ngettext(sum(huh), "is", "are"), "not defined"), call.=FALSE)
  df <- SparseEntries(X)
  # discard other indices
  nonmargin <- setdiff(seq_along(dimX), MARGIN)
  df <- df[ , -nonmargin, drop=FALSE]
  # implicitly accumulate
  result <- EntriesToSparse(df, dimX[MARGIN])
  return(result)
}

tensor1x1 <- function(A, B) {
  ## equivalent of tensor(A, B, 1, 1)
  ## when A is a vector and B is a sparse array.
  stopifnot(length(dim(B)) == 3)
  A <- as.vector(as.matrix(A))
  stopifnot(length(A) == dim(B)[1])
  if(is.array(B)) {
    result <- tensor::tensor(A,B,1,1)
  } else if(inherits(B, "sparse3Darray")) {
    result <- sparseMatrix(i=B$j,
                           j=B$k,
                           x=B$x * A[B$i], # values for same (i,j) are summed
                           dims=dim(B)[-1],
                           dimnames=dimnames(B)[2:3])
  } else stop("Format of B not understood", call.=FALSE)
  return(result)
}

tenseur <- local({

  tenseur <- function(A, B, alongA=integer(0), alongB=integer(0)) {
    #' full arrays?
    if(isfull(A) && isfull(B))
      return(tensor::tensor(A=A, B=B, alongA=alongA, alongB=alongB))
    #' check dimensions
    dimA <- dim(A) %orifnull% length(A)
    dnA <- dimnames(A)
    if(is.null(dnA))
      dnA <- rep(list(NULL), length(dimA))
    dimB <- dim(B) %orifnull% length(B)
    dnB <- dimnames(B)
    if(is.null(dnB))
      dnB <- rep(list(NULL), length(dimB))
    #' check 'along'
    if (length(alongA) != length(alongB)) 
      stop("\"along\" vectors must be same length")
    mtch <- dimA[alongA] == dimB[alongB]
    if (any(is.na(mtch)) || !all(mtch)) 
      stop("Mismatch in \"along\" dimensions")
    #' dimensions of result
    retainA <- !(seq_along(dimA) %in% alongA)
    retainB <- !(seq_along(dimB) %in% alongB)
    dimC <- c(dimA[retainA], dimB[retainB])
    nC <- length(dimC)
    if(nC > 3)
      stop("Sorry, sparse arrays of more than 3 dimensions are not supported",
           call.=FALSE)
    #' extract indices and values of nonzero entries
    dfA <- SparseEntries(A)
    dfB <- SparseEntries(B)
    #' assemble all tuples which contribute 
    if(length(alongA) == 0) {
      #' outer product
      dfC <- outersparse(dfA, dfB)
    } else {
      if(length(alongA) == 1) {
        Acode <- dfA[,alongA]
        Bcode <- dfB[,alongB]
      } else {
        Along <- unname(as.list(dfA[,alongA, drop=FALSE]))
        Blong <- unname(as.list(dfB[,alongB, drop=FALSE]))
        Acode <- do.call(paste, append(Along, list(sep=",")))
        Bcode <- do.call(paste, append(Blong, list(sep=",")))
      }
      lev <- unique(c(Acode,Bcode))
      Acode <- factor(Acode, levels=lev)
      Bcode <- factor(Bcode, levels=lev)
      splitA <- split(dfA, Acode)
      splitB <- split(dfB, Bcode)
      splitC <- mapply(outersparse, splitA, splitB, SIMPLIFY=FALSE)
      dfC <- rbindCompatibleDataFrames(splitC)
    }
    #' form product of contributing entries
    dfC$x <- with(dfC, A.x * B.x)
    #' retain only appropriate columns
    retain <- c(retainA, FALSE, retainB, FALSE, TRUE)
    dfC <- dfC[, retain, drop=FALSE]
    #' collect result
    result <- EntriesToSparse(dfC, dimC)
    return(result)
  }

  isfull <- function(z) {
    if(is.array(z) || is.matrix(z) || is.data.frame(z)) return(TRUE)
    if(inherits(z, c("sparseVector", "sparseMatrix", "sparse3Darray")))
      return(FALSE)
    return(TRUE)
  }
  
  outersparse <- function(dfA, dfB) {
    if(is.null(dfA) || is.null(dfB)) return(NULL)
    IJ <- expand.grid(I=seq_len(nrow(dfA)),
                      J=seq_len(nrow(dfB)))
    dfC <- with(IJ, cbind(A=dfA[I,,drop=FALSE], B=dfB[J,,drop=FALSE]))
    return(dfC)
  }

  tenseur
})

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
            y  = as.double(numeric(m * m)),
            PACKAGE = "spatstat")
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
            y  = as.double(numeric(m * m)),
            PACKAGE = "spatstat")
  }
  y <- matrix(z$y, m, m)
  dimnames(y) <- rep(dimnames(x)[1], 2)
  return(y)
}

