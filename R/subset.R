##
## subset.R
##
## Methods for 'subset'
##
##   $Revision: 1.6 $  $Date: 2018/10/02 03:50:44 $

subset.ppp <- function(x, subset, select, drop=FALSE, ...) {
  stopifnot(is.ppp(x))
  w <- as.owin(x)
  y <- as.data.frame(x)
  r <- if (missing(subset)) {
    rep_len(TRUE, nrow(y))
  } else {
    e <- substitute(subset)
    r <- eval(e, y, parent.frame())
    if(!is.logical(r))
      r <- ppsubset(x, r, "subset", fatal=TRUE)
    r & !is.na(r)
  }
  vars <- if (missing(select)) {
    TRUE
  } else {
    ## create an environment in which column names are mapped to their positions
    nl <- as.list(seq_along(y))
    names(nl) <- names(y)
    if(length(nl) > 3) {
      ## multiple columns of marks: add the name 'marks'
      nl <- append(nl, list(marks=3:length(nl)))
    }
    eval(substitute(select), nl, parent.frame())
  }
  ## ensure columns include coordinates
  nama <- names(y)
  names(nama) <- nama
  vars <- union(c("x", "y"), nama[vars])
  ## take subset
  z <- y[r, vars, drop = FALSE]
  ## reinstate as point pattern
  out <- as.ppp(z, W=w, check=FALSE)
  if(drop)
    out <- out[drop=TRUE]
  return(out)
}

subset.pp3 <- subset.lpp <- subset.ppx <- function(x, subset, select, drop=FALSE, ...) {
  y <- as.data.frame(x)
  r <- if (missing(subset)) 
    rep_len(TRUE, nrow(y))
  else {
    e <- substitute(subset)
    r <- eval(e, y, parent.frame())
    if(!is.logical(r))
      r <- ppsubset(x, r, "subset", fatal=TRUE)
    r & !is.na(r)
  }
  vars <- if (missing(select)) 
    TRUE
  else {
    ## create an environment in which column names are mapped to their positions
    nl <- as.list(seq_along(y))
    names(nl) <- names(y)
    if(!("marks" %in% names(y)) && any(ismark <- (x$ctype == "mark"))) {
      ## add the symbol 'marks' 
      nl <- append(nl, list(marks=which(ismark)))
    }
    eval(substitute(select), nl, parent.frame())
  }
  ## ensure columns include coordinates 
  nama <- names(y)
  names(nama) <- nama
  vars <- union(names(coords(x)), nama[vars])
  ## take subset
  z <- y[r, vars, drop = FALSE]
  ## reinstate as point pattern
  ctype <- as.character(x$ctype)[match(vars, nama)]
  out <- ppx(z, domain=x$domain, coord.type=ctype)
  ## drop unused factor levels
  if(drop)
    out <- out[drop=TRUE]
  ## reinstate class
  class(out) <- class(x)
  return(out)
}

subset.psp <- function(x, subset, select, drop=FALSE, ...) {
  stopifnot(is.psp(x))
  w <- Window(x)
  y <- as.data.frame(x)
  r <- if (missing(subset)) {
    rep_len(TRUE, nrow(y))
  } else {
    e <- substitute(subset)
    r <- eval(e, y, parent.frame())
    if(!is.logical(r))
      stop("Argument 'subset' should be a logical vector", call.=FALSE)
    r & !is.na(r)
  }
  vars <- if (missing(select)) {
    TRUE
  } else {
    ## create an environment in which column names are mapped to their positions
    nl <- as.list(seq_along(y))
    names(nl) <- names(y)
    if(length(nl) > 3) {
      ## multiple columns of marks: add the name 'marks'
      nl <- append(nl, list(marks=3:length(nl)))
    }
    eval(substitute(select), nl, parent.frame())
  }
  ## ensure columns include coordinates
  nama <- names(y)
  names(nama) <- nama
  vars <- union(c("x0", "y0", "x1", "y1"), nama[vars])
  ## take subset
  z <- y[r, vars, drop = FALSE]
  ## reinstate as line segment pattern
  out <- as.psp(z, window=w, check=FALSE)
  if(drop)
    out <- out[drop=TRUE]
  return(out)
}
