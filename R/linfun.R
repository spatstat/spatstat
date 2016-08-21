#
#   linfun.R
#
#   Class of functions of location on a linear network
#
#   $Revision: 1.10 $   $Date: 2016/08/21 03:52:39 $
#

linfun <- function(f, L) {
  stopifnot(is.function(f))
  stopifnot(inherits(L, "linnet"))
  fargs <- names(formals(f))
  needargs <- c("x", "y", "seg", "tp")
  if(!all(needargs %in% fargs))
    stop(paste("Function must have formal arguments",
               commasep(sQuote(needargs))),
         call.=FALSE)
  otherfargs <- setdiff(fargs, needargs)
  g <- function(...) {
    argh <- list(...)
    extra <- names(argh) %in% otherfargs
    if(!any(extra)) {
      X <- as.lpp(..., L=L)
      value <- do.call(f, as.list(coords(X)))
    } else {
      extrargs <- argh[extra]
      mainargs <- argh[!extra]
      X <- do.call(as.lpp, append(mainargs, list(L=L)))
      value <- do.call(f, append(as.list(coords(X)), extrargs))
    }
    return(value)
  }
  class(g) <- c("linfun", class(g))
  attr(g, "L") <- L
  attr(g, "f") <- f
  return(g)
}

print.linfun <- function(x, ...) {
  L <- as.linnet(x)
  if(!is.null(explain <- attr(x, "explain"))) {
    explain(x)
  } else {
    splat("Function on linear network:")
    print(attr(x, "f"), ...)
    splat("Function domain:")
    print(L)
  }
  invisible(NULL)
}

summary.linfun <- function(object, ...) { print(object, ...) }

as.linim.linfun <- function(X, L, ..., eps = NULL, dimyx = NULL, xy = NULL) {
  if(missing(L) || is.null(L))
    L <- as.linnet(X)
  # create template
  Y <- as.linim(1, L, eps=eps, dimyx=dimyx, xy=xy)
  # extract coordinates of sample points along network
  df <- attr(Y, "df")
  coo <- df[, c("x", "y", "mapXY", "tp")]
  colnames(coo)[3] <- "seg"
  # evaluate function at sample points
  vals <- do.call(X, append(as.list(coo), list(...)))
  # write values in data frame
  df$values <- vals
  # write values in pixel array
  typ <- typeof(vals)
  storage.mode(Y$v) <- typ
  Y[] <- vals
  Y$type <- if(typ == "double") "real" else typ
  attr(Y, "df") <- df
  return(Y)
}

as.data.frame.linfun <- function(x, ...) {
  as.data.frame(as.linim(x, ...))
}

as.linfun.linim <- function(X, ...) {
  trap.extra.arguments(..., .Context="as.linfun.linim")
  ## extract info
  L <- as.linnet(X)
  df <- attr(X, "df")
  ## function values and corresponding locations
  values <- df$values
  locations <- with(df, as.lpp(x=x, y=y, seg=mapXY, tp=tp, L=L))
  ## Function that maps any spatial location to the nearest data location 
  nearestloc <- nnfun(locations)
  ## Function that reads value at nearest data location
  f <- function(x, y, seg, tp) {
    values[nearestloc(x,y,seg,tp)]
  }
  g <- linfun(f, L)
  return(g)
}

plot.linfun <- function(x, ..., L=NULL, eps = NULL, dimyx = NULL, xy = NULL,
                        main) {
  if(missing(main)) main <- short.deparse(substitute(x))
  if(is.null(L)) L <- as.linnet(x)
  argh <- list(...)
  otherfargs <- get("otherfargs", envir=environment(x))
  xtra <- names(argh) %in% otherfargs
  if(!any(xtra)) {
    Z <- as.linim(x, eps=eps, dimyx=dimyx, xy=xy, L=L)
    rslt <- plot(Z, ..., main=main)
  } else {
    Z <- do.call(as.linim, append(list(x, eps=eps, dimyx=dimyx, xy=xy, L=L),
                                  argh[xtra]))
    rslt <- do.call(plot.linim, append(list(Z, main=main), argh[!xtra]))
  }
  return(invisible(rslt))
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

integral.linfun <- function(f, domain=NULL, ...) {
  integral(as.linim(f), domain=domain, ...)
}

as.linfun <- function(X, ...) {
  UseMethod("as.linfun")
}

