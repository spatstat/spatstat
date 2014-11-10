#
#   linfun.R
#
#   Class of functions of location on a linear network
#
#   $Revision: 1.1 $   $Date: 2014/10/24 00:22:30 $
#

linfun <- function(f, L) {
  stopifnot(is.function(f))
  stopifnot(inherits(L, "linnet"))
  needargs <- c("x", "y", "seg", "tp")
  if(!all(needargs %in% names(formals(f))))
    stop(paste("f must have arguments named", commasep(sQuote(needargs))))
  class(f) <- c("linfun", class(f))
  attr(f, "L") <- L
  return(f)
}

print.linfun <- function(x, ...) {
  L <- as.linnet(x)
  if(!is.null(explain <- attr(x, "explain"))) {
    explain(x)
  } else {
    cat("Function on linear network\n")
    print(as.function(x), ...)
    cat("Function domain:\n")
    print(L)
  }
  invisible(NULL)
}

as.linim.linfun <- function(X, L, ..., eps = NULL, dimyx = NULL, xy = NULL) {
  if(missing(L) || is.null(L))
    L <- as.linnet(X)
  # create template
  Y <- as.linim(1, L, eps=eps, dimyx=dimyx, xy=xy)
  # extract (x,y) and local coordinates
  df <- attr(Y, "df")
  coo <- df[, c("x", "y", "mapXY", "tp")]
  colnames(coo)[3] <- "seg"
  # evaluate function
  vals <- do.call(X, append(as.list(coo), list(...)))
  # replace values
  df$values <- vals
  attr(Y, "df") <- df
  Y[!is.na(Y$v)] <- vals
  return(Y)
}
  
plot.linfun <- function(x, ..., L=NULL, eps = NULL, dimyx = NULL, xy = NULL,
                        main) {
  if(missing(main)) main <- short.deparse(substitute(x))
  if(is.null(L)) L <- as.linnet(x)
  Z <- as.linim(x, eps=eps, dimyx=dimyx, xy=xy, L=L)
  plot(Z, ..., main=main)
}

as.owin.linfun <- function(W, ...) {
  as.owin(as.linnet(W))
}

as.linnet.linfun <- function(X, ...) {
  attr(X, "L")
}

as.function.linfun <- function(x, ...) {
  nax <- names(attributes(x))
  if(!is.null(nax)) {
    retain <- (nax == "srcref")
    attributes(x)[!retain] <- NULL
  }
  return(x)
}

