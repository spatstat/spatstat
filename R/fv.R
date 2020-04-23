##
##
##    fv.R
##
##    class "fv" of function value objects
##
##    $Revision: 1.159 $   $Date: 2020/04/23 04:04:41 $
##
##
##    An "fv" object represents one or more related functions
##    of the same argument, such as different estimates of the K function.
##
##    It is a data.frame with additional attributes
##    
##         argu       column name of the function argument (typically "r")
##
##         valu       column name of the recommended function
##
##         ylab       generic label for y axis e.g. K(r)
##
##         fmla       default plot formula
##
##         alim       recommended range of function argument
##
##         labl       recommended xlab/ylab for each column
##
##         desc       longer description for each column 
##
##         unitname   name of unit of length for 'r'
##
##         shade      (optional) column names of upper & lower limits
##                    of shading - typically a confidence interval
##
##    Objects of this class are returned by Kest(), etc
##
##################################################################
## creator

fv <- function(x, argu="r", ylab=NULL, valu, fmla=NULL,
               alim=NULL, labl=names(x), desc=NULL, unitname=NULL,
               fname=NULL, yexp=ylab) {
  stopifnot(is.data.frame(x))
  ## check arguments
  stopifnot(is.character(argu))
  if(!is.null(ylab))
    stopifnot(is.character(ylab) || is.language(ylab))
  if(!missing(yexp)) {
    if(is.null(yexp)) yexp <- ylab
    else stopifnot(is.language(yexp))
  }
  stopifnot(is.character(valu))
  
  if(!(argu %in% names(x)))
    stop(paste(sQuote("argu"), "must be the name of a column of x"))

  if(!(valu %in% names(x)))
    stop(paste(sQuote("valu"), "must be the name of a column of x"))

  if(is.null(fmla))
    fmla <- paste(valu, "~", argu)
  else if(inherits(fmla, "formula")) {
    ## convert formula to string
    fmla <- flat.deparse(fmla)
  } else if(!is.character(fmla))
    stop(paste(sQuote("fmla"), "should be a formula or a string"))

  if(missing(alim)) {
    ## Note: if alim is given as NULL, it is not changed.
    argue <- x[[argu]]
    alim <- range(argue[is.finite(argue)])
  } else if(!is.null(alim)) {
    if(!is.numeric(alim) || length(alim) != 2)
      stop(paste(sQuote("alim"), "should be a vector of length 2"))
  }
  
  if(!is.character(labl))
    stop(paste(sQuote("labl"), "should be a vector of strings"))
  stopifnot(length(labl) == ncol(x))
  if(is.null(desc))
    desc <- character(ncol(x))
  else {
    stopifnot(is.character(desc))
    stopifnot(length(desc) == ncol(x))
    nbg <- is.na(desc)
    if(any(nbg)) desc[nbg] <- ""
  }
  if(!is.null(fname))
    stopifnot(is.character(fname) && length(fname) %in% 1:2)
  ## pack attributes
  attr(x, "argu") <- argu
  attr(x, "valu") <- valu
  attr(x, "ylab") <- ylab
  attr(x, "yexp") <- yexp
  attr(x, "fmla") <- fmla
  attr(x, "alim") <- alim
  attr(x, "labl") <- labl
  attr(x, "desc") <- desc
  attr(x, "units") <- as.unitname(unitname)
  attr(x, "fname") <- fname
  attr(x, "dotnames") <- NULL
  attr(x, "shade") <- NULL
  ## 
  class(x) <- c("fv", class(x))
  return(x)
}

.Spatstat.FvAttrib <- c(
                        "argu",
                        "valu",
                        "ylab",
                        "yexp",
                        "fmla",
                        "alim",
                        "labl",
                        "desc",
                        "units",
                        "fname",
                        "dotnames",
                        "shade")

as.data.frame.fv <- function(x, ...) {
  stopifnot(is.fv(x))
  fva <- .Spatstat.FvAttrib
  attributes(x)[fva] <- NULL
  class(x) <- "data.frame"
  x
}

is.fv <- function(x) {
  inherits(x, "fv")
}

## 

as.fv <- function(x) { UseMethod("as.fv") }

as.fv.fv <- function(x) x

as.fv.data.frame <- function(x) {
  if(ncol(x) < 2) stop("Need at least 2 columns")
  return(fv(x, names(x)[1L], , names(x)[2L]))
}

as.fv.matrix <- function(x)  {
  y <- as.data.frame(x)
  if(any(bad <- is.na(names(y))))
    names(y)[bad] <- paste0("V", which(bad))
  return(as.fv.data.frame(y))
}

## other methods for as.fv are described in the files for the relevant classes.

vanilla.fv <- function(x) {
  ## remove everything except basic fv characteristics
  retain <- c("names", "row.names", .Spatstat.FvAttrib)
  attributes(x) <- attributes(x)[retain]
  class(x) <- c("fv", "data.frame")
  return(x)
}

print.fv <- local({
  
  maxwords <- function(z, m) { max(0, which(cumsum(nchar(z) + 1) <= m+1)) }
  usewords <- function(z, n) paste(z[1:n], collapse=" ")

  print.fv <- function(x, ..., tight=FALSE) {
    verifyclass(x, "fv")
    terselevel <- spatstat.options("terse")
    showlabels <- waxlyrical('space', terselevel)
    showextras <- waxlyrical('extras', terselevel)
    nama <- names(x)
    a <- attributes(x)
    if(!is.null(ylab <- a$ylab)) {
      if(is.language(ylab))
        ylab <- flat.deparse(ylab)
    }
    if(!inherits(x, "envelope")) {
      splat("Function value object",
            paren(paste("class", sQuote("fv"))))
      if(!is.null(ylab)) {
        xlab <- fvlabels(x, expand=TRUE)[[a$argu]]
        splat("for the function", xlab, "->", ylab)
      }
    }
    ## Descriptions ..
    desc <- a$desc
    ## .. may require insertion of ylab
    if(!is.null(ylab))
      desc <- sprintf(desc, ylab)
    ## Labels ..
    labl <- fvlabels(x, expand=TRUE)
    ## Avoid overrunning text margin
    maxlinewidth <- options('width')[[1L]]
    key.width <- max(nchar(nama))
    labl.width <- if(showlabels) max(nchar(labl), nchar("Math.label")) else 0
    desc.width <- max(nchar(desc), nchar("Description"))
    fullwidth <- key.width + labl.width + desc.width + 2
    if(fullwidth > maxlinewidth && tight) {
      ## try shortening the descriptions so that it all fits on one line
      spaceleft <- maxlinewidth - (key.width + labl.width + 2)
      desc <- truncline(desc, spaceleft)
      desc.width <- max(nchar(desc), nchar("Description"))    
      fullwidth <- key.width + labl.width + desc.width + 2
    }
    spaceleft <- maxlinewidth - (key.width + 1)
    if(desc.width > spaceleft) {
      ## Descriptions need to be truncated to max line width
      desc <- truncline(desc, spaceleft)
      desc.width <- max(nchar(desc), nchar("Description"))    
      fullwidth <- key.width + labl.width + desc.width + 2
    }
    if(showextras) {
      fullwidth <- pmin(maxlinewidth, fullwidth)
      fullline <- paste0(rep(".", fullwidth), collapse="")
      cat(fullline, fill=TRUE)
    }
    df <- data.frame(Math.label=labl,
                     Description=desc,
                     row.names=nama,
                     stringsAsFactors=FALSE)
    if(!showlabels) df <- df[,-1L,drop=FALSE]
    print(df, right=FALSE)
  ##
    if(showextras) {
      cat(fullline, fill=TRUE)
      splat("Default plot formula: ",
            flat.deparse(as.formula(a$fmla)))
      splat("where", dQuote("."), "stands for",
            commasep(sQuote(fvnames(x, ".")), ", "))
      if(length(a$shade)) 
        splat("Columns", commasep(sQuote(a$shade)), 
              "will be plotted as shading (by default)")
      alim <- a$alim
      splat("Recommended range of argument",
            paste0(a$argu, ":"),
            if(!is.null(alim)) prange(signif(alim, 5)) else "not specified")
      rang <- signif(range(with(x, .x)), 5)
      splat("Available range", "of argument",
            paste0(a$argu, ":"), prange(rang))
      ledge <- summary(unitname(x))$legend
      if(!is.null(ledge))
        splat(ledge)
    }
    return(invisible(NULL))
  }

  print.fv
})



## manipulating the names in fv objects

.Spatstat.FvAbbrev <- c(
                        ".x",
                        ".y",
                        ".s",
                        ".",
                        "*",
                        ".a")

fvnames <- function(X, a=".") {
  verifyclass(X, "fv")
  if(!is.character(a))
    stop("argument a must be a character string")
  if(length(a) != 1) return(lapply(a, function(b, Z) fvnames(Z, b), Z=X))
  namesX <- names(X)
  if(a %in% namesX) return(a)
  vnames <- setdiff(namesX, attr(X, "argu"))
  answer <- switch(a,
                   ".y" = attr(X, "valu"),
                   ".x" = attr(X, "argu"),
                   ".s" = attr(X, "shade"),
                   ".a" = vnames,
                   "*"  = rev(vnames),
                   "."  = attr(X, "dotnames") %orifnull% rev(vnames),
                   {
                     stop(paste("Unrecognised abbreviation", sQuote(a)),
                          call.=FALSE)
                     })
  return(answer)
}

"fvnames<-" <- function(X, a=".", value) {
  verifyclass(X, "fv")
  if(!is.character(a) || length(a) > 1)
    stop(paste("argument", sQuote("a"), "must be a character string"))
  ## special cases
  if(a == "." && length(value) == 0) {
    ## clear the dotnames
    attr(X, "dotnames") <- NULL
    return(X)
  }
  if(a == ".a" || a == "*") {
    warning("Nothing changed; use names(X) <- value to change names",
            call.=FALSE)
    return(X)
  }
  ## validate the names
  switch(a,
         ".x"=,
         ".y"={
           if(!is.character(value) || length(value) != 1)
             stop("value should be a single string")
         },
         ".s"={
           if(!is.character(value) || length(value) != 2)
             stop("value should be a vector of 2 character strings")
         },
         "."={
           if(!is.character(value))
             stop("value should be a character vector")
         },
         stop(paste("Unrecognised abbreviation", dQuote(a)))
       )
  ## check the names match existing column names
  tags <- names(X)
  if(any(nbg <- !(value %in% tags))) 
    stop(paste(ngettext(sum(nbg), "The string", "The strings"),
               commasep(dQuote(value[nbg])),
               ngettext(sum(nbg),
                        "does not match the name of any column of X", 
                        "do not match the names of any columns of X")))
  ## reassign names
  switch(a,
         ".x"={
           attr(X, "argu") <- value
         },
         ".y"={
           attr(X, "valu") <- value
         },
         ".s"={
           attr(X, "shade") <- value
         },
         "."={
           attr(X, "dotnames") <- value
         })
  return(X)
}

"names<-.fv" <- function(x, value) {
  nama <- colnames(x)
  indx <- which(nama == fvnames(x, ".x"))
  indy <- which(nama == fvnames(x, ".y"))
  inds <- which(nama %in% fvnames(x, ".s"))
  ind. <- which(nama %in% fvnames(x, "."))
  ## rename columns of data frame
  x <- NextMethod("names<-")
  ## adjust other tags
  fvnames(x, ".x") <- value[indx]
  fvnames(x, ".y") <- value[indy]
  fvnames(x, ".")  <- value[ind.]
  if(length(inds))
    fvnames(x, ".s") <- value[inds]
  namemap <- setNames(lapply(value, as.name), nama)
  formula(x) <- flat.deparse(eval(substitute(substitute(fom, um),
                                             list(fom=as.formula(formula(x)),
                                                  um=namemap))))
  return(x)
}

fvlabels <- function(x, expand=FALSE) {
  lab <- attr(x, "labl")
  if(expand && !is.null(fname <- attr(x, "fname"))) {
    ## expand plot labels using function name
    nstrings <- max(substringcount("%s", lab))
    ## pad with blanks
    nextra <- nstrings - length(fname)
    if(nextra > 0) 
      fname <- c(fname, rep("", nextra))
    ## render
    lab <- do.call(sprintf, append(list(lab), as.list(fname)))
  }
  ## remove empty space
  lab <- gsub(" ", "", lab)
  names(lab) <- names(x)
  return(lab)
}

"fvlabels<-" <- function(x, value) {
  stopifnot(is.fv(x))
  stopifnot(is.character(value))
  stopifnot(length(value) == length(fvlabels(x)))
  attr(x, "labl") <- value
  return(x)
}

flatfname <- function(x) {
  fn <- if(is.character(x)) x else attr(x, "fname")
  if(length(fn) > 1)
    fn <- paste0(fn[1L], "[", paste(fn[-1L], collapse=" "), "]")
  as.name(fn)
}

makefvlabel <- function(op=NULL, accent=NULL, fname, sub=NULL, argname="r") {
  ## de facto standardised label
  a <- "%s"
  if(!is.null(accent)) 
    a <- paste0(accent, paren(a))     ## eg hat(%s)
  if(!is.null(op))
    a <- paste0("bold", paren(op), "~", a)  ## eg bold(var)~hat(%s)
  if(is.null(sub)) {
    if(length(fname) != 1) {
      a <- paste0(a, "[%s]")
      a <- paren(a, "{")
    }
  } else {
    if(length(fname) == 1) {
      a <- paste0(a, paren(sub, "["))
    } else {
      a <- paste0(a, paren("%s", "["), "^", paren(sub, "{"))
      a <- paren(a, "{")
    }
  } 
  a <- paste0(a, paren(argname))
  return(a)
}

fvlabelmap <- local({
  magic <- function(x) {
    subx <- paste("substitute(", x, ", NULL)")
    out <- try(eval(parse(text=subx)), silent=TRUE)
    if(inherits(out, "try-error"))
      out <- as.name(make.names(subx))
    out
  }

  fvlabelmap <- function(x, dot=TRUE) {
    labl <- fvlabels(x, expand=TRUE)
    ## construct mapping from identifiers to labels
    map <- as.list(labl)
    map <- lapply(map, magic)
    names(map) <- colnames(x)
    if(dot) {
      ## also map "." and ".a" to name of target function
      if(!is.null(ye <- attr(x, "yexp")))
        map <- append(map, list("."=ye, ".a"=ye))
      ## map other fvnames to their corresponding labels
      map <- append(map, list(".x"=map[[fvnames(x, ".x")]],
                              ".y"=map[[fvnames(x, ".y")]]))
      if(length(fvnames(x, ".s"))) {
        shex <- unname(map[fvnames(x, ".s")])
        shadexpr <- substitute(c(A,B), list(A=shex[[1L]], B=shex[[2L]]))
        map <- append(map, list(".s" = shadexpr))
      }
    }
    return(map)
  }

  fvlabelmap
})

## map from abbreviations to expressions involving the column names,
## for use in eval(substitute(...))
fvexprmap <- function(x) {
  dotnames <- fvnames(x, ".")
  u <- if(length(dotnames) == 1) as.name(dotnames) else 
       as.call(lapply(c("cbind", dotnames), as.name))
  ux <- as.name(fvnames(x, ".x"))
  uy <- as.name(fvnames(x, ".y"))
  umap <- list(.=u, .a=u, .x=ux, .y=uy)
  if(length(shnm <- fvnames(x, ".s"))) {
    shadexpr <- substitute(cbind(A,B), list(A=as.name(shnm[1L]),
                                            B=as.name(shnm[2L])))
    umap <- append(umap, list(.s = shadexpr))
  }
  return(umap)
}

fvlegend <- local({

  fvlegend <- function(object, elang) {
    ## Compute mathematical legend(s) for column(s) in fv object 
    ## transformed by language expression 'elang'.
    ## The expression must already be in 'expanded' form.
    ## The result is an expression, or expression vector.
    ## The j-th entry of the vector is an expression for the
    ## j-th column of function values.
    ee <- distributecbind(as.expression(elang))
    map <- fvlabelmap(object, dot = TRUE)
    eout <- as.expression(lapply(ee, invokemap, map=map))
    return(eout)
  }

  invokemap <- function(ei, map) {
    eval(substitute(substitute(e, mp), list(e = ei, mp = map)))
  }
  
  fvlegend
})


bind.fv <- function(x, y, labl=NULL, desc=NULL, preferred=NULL, clip=FALSE) {
  verifyclass(x, "fv")
  ax <- attributes(x)
  if(is.fv(y)) {
    ## y is already an fv object
    ay <- attributes(y)
    if(!identical(ax$fname, ay$fname)) {
      ## x and y represent different functions
      ## expand the labels separately 
      fvlabels(x) <- fvlabels(x, expand=TRUE)
      fvlabels(y) <- fvlabels(y, expand=TRUE)
      ax <- attributes(x)
      ay <- attributes(y)
    }
    ## check compatibility of 'r' values
    xr <- ax$argu
    yr <- ay$argu
    rx <- x[[xr]]
    ry <- y[[yr]]
    if(length(rx) != length(ry)) {
      if(!clip) 
        stop("fv objects x and y have incompatible domains")
      # restrict both objects to a common domain
      ra <- intersect.ranges(range(rx), range(ry))
      x <- x[inside.range(rx, ra), ]
      y <- y[inside.range(ry, ra), ]
      rx <- x[[xr]]
      ry <- y[[yr]]
    }
    if(length(rx) != length(ry) || max(abs(rx-ry)) > .Machine$double.eps)
      stop("fv objects x and y have incompatible values of r")
    ## reduce y to data frame and strip off 'r' values
    ystrip <- as.data.frame(y)
    yrpos <- which(colnames(ystrip) == yr)
    ystrip <- ystrip[, -yrpos, drop=FALSE]
    ## determine descriptors
    if(is.null(labl)) labl <- attr(y, "labl")[-yrpos]
    if(is.null(desc)) desc <- attr(y, "desc")[-yrpos]
    ##
    y <- ystrip
  } else {
    ## y is a matrix or data frame
    y <- as.data.frame(y)
  }
  
  ## check for duplicated column names
  allnames <- c(colnames(x), colnames(y))
  if(any(dup <- duplicated(allnames))) {
    nbg <- unique(allnames[dup])
    nn <- length(nbg)
    warning(paste("The column",
                  ngettext(nn, "name", "names"),
                  commasep(sQuote(nbg)),
                  ngettext(nn, "was", "were"),
                  "duplicated. Unique names were generated"))
    allnames <- make.names(allnames, unique=TRUE, allow_ = FALSE)
    colnames(y) <- allnames[ncol(x) + seq_len(ncol(y))]
  }
      
  if(is.null(labl))
    labl <- paste("%s[", colnames(y), "](r)", sep="")
  else if(length(labl) != ncol(y))
    stop(paste("length of", sQuote("labl"),
               "does not match number of columns of y"))
  if(is.null(desc))
    desc <- character(ncol(y))
  else if(length(desc) != ncol(y))
    stop(paste("length of", sQuote("desc"),
               "does not match number of columns of y"))
  if(is.null(preferred))
    preferred <- ax$valu

  xy <- cbind(as.data.frame(x), y)
  z <- fv(xy, ax$argu, ax$ylab, preferred, ax$fmla, ax$alim,
          c(ax$labl, labl),
          c(ax$desc, desc),
          unitname=unitname(x),
          fname=ax$fname,
          yexp=ax$yexp)
  return(z)
}

cbind.fv <- function(...) {
  a <- list(...)
  n <- length(a)
  if(n == 0)
    return(NULL)
  if(n == 1) {
    ## single argument - extract it
    a <- a[[1L]]
    ## could be an fv object 
    if(is.fv(a))
      return(a)
    n <- length(a)
  }
  z <- a[[1L]]
  if(!is.fv(z))
    stop("First argument should be an object of class fv")
  if(n > 1)
    for(i in 2:n) 
      z <- bind.fv(z, a[[i]])
  return(z)
}

collapse.anylist <-
collapse.fv <- local({

  collapse.fv <- function(object, ..., same=NULL, different=NULL) {
    if(is.fv(object)) {
      x <- list(object, ...)
    } else if(inherits(object, "anylist")) {
      x <- append(object, list(...))
    } else if(is.list(object) && all(sapply(object, is.fv))) {
      x <- append(object, list(...))
    } else stop("Format not understood")
    if(!all(unlist(lapply(x, is.fv))))
      stop("arguments should be objects of class fv")
    if(is.null(same)) same <- character(0)
    if(is.null(different)) different <- character(0)
    if(anyDuplicated(c(same, different)))
      stop(paste("The arguments", sQuote("same"), "and", sQuote("different"),
                 "should not have entries in common"))
    either <- c(same, different)
    ## validate
    if(length(either) == 0)
      stop(paste("At least one column of values must be selected",
                 "using the arguments", sQuote("same"), "and",
                 sQuote("different")))
    nbg <- unique(unlist(lapply(x, missingnames, expected=either)))
    if((nbad <- length(nbg)) > 0)
      stop(paste(ngettext(nbad, "The name", "The names"),
                 commasep(sQuote(nbg)),
                 ngettext(nbad, "is", "are"),
                 "not present in the function objects"))
    ## names for different versions
    versionnames <- names(x)
    if(is.null(versionnames))
      versionnames <- paste("x", seq_along(x), sep="")
    shortnames <- abbreviate(versionnames, minlength=12)
    ## extract the common values
    y <- x[[1L]]
    xname <- fvnames(y, ".x")
    yname <- fvnames(y, ".y")
    if(length(same) == 0) {
      ## The column of 'preferred values' .y cannot be deleted
      ## retain .y for now and delete it later.
      z <- y[, c(xname, yname)]
    } else {
      same <-  fvnames(y, same) # expand abbreviations if present
      if(!(yname %in% same))
        fvnames(y, ".y") <- same[1L]
      z <- y[, c(xname, same)]
    }
    dotnames <- same
    ## now merge the different values
    if(length(different)) {
      for(i in seq_along(x)) {
        ## extract values for i-th object
        xi <- x[[i]]
        diffi <- fvnames(xi, different) # expand abbreviations if present
        wanted <- (names(xi) %in% diffi)
        if(any(wanted)) {
          y <- as.data.frame(xi)[, wanted, drop=FALSE]
          desc <- attr(xi, "desc")[wanted]
          labl <- attr(xi, "labl")[wanted]
          ## relabel
          prefix <- shortnames[i]
          preamble <- versionnames[i]
          names(y) <- if(ncol(y) == 1) prefix else paste(prefix,names(y),sep="")
          dotnames <- c(dotnames, names(y))
          ## glue onto fv object
          z <- bind.fv(z, y,
                       labl=paste(prefix, labl, sep="~"),
                       desc=paste(preamble, desc))
        }
      }
    }
    if(length(same) == 0) {
      ## remove the second column which was retained earlier
      fvnames(z, ".y") <- names(z)[3L]
      z <- z[, -2L]
    }
    fvnames(z, ".") <- dotnames
    return(z)
  }

  
  missingnames <- function(z, expected) {
    known <- c(names(z), .Spatstat.FvAbbrev)
    absent <- is.na(match(expected, known)) 
    return(expected[absent])
  }
  
  collapse.fv
})

## rename one of the columns of an fv object
tweak.fv.entry <- function(x, current.tag, new.labl=NULL, new.desc=NULL, new.tag=NULL) {
  hit <- (names(x) == current.tag)
  if(!any(hit))
    return(x)
  ## update descriptions of column
  i <- min(which(hit))
  if(!is.null(new.labl)) attr(x, "labl")[i] <- new.labl
  if(!is.null(new.desc)) attr(x, "desc")[i] <- new.desc
  ## adjust column tag
  if(!is.null(new.tag)) {
    names(x)[i] <- new.tag
    ## update dotnames
    dn <- fvnames(x, ".")
    if(current.tag %in% dn ) {
      dn[dn == current.tag] <- new.tag
      fvnames(x, ".") <- dn
    }
    ## if the tweaked column is the preferred value, adjust accordingly
    if(attr(x, "valu") == current.tag)
      attr(x, "valu") <- new.tag
    ## if the tweaked column is the function argument, adjust accordingly
    if(attr(x, "argu") == current.tag)
      attr(x, "valu") <- new.tag
  }
  return(x)
}


## change some or all of the auxiliary text in an fv object
rebadge.fv <- function(x, new.ylab, new.fname,
                       tags, new.desc, new.labl,
                       new.yexp=new.ylab, new.dotnames,
                       new.preferred, new.formula, new.tags) {
  if(!missing(new.ylab)) 
    attr(x, "ylab") <- new.ylab
  if(!missing(new.yexp) || !missing(new.ylab))
    attr(x, "yexp") <- new.yexp
  if(!missing(new.fname))
    attr(x, "fname") <- new.fname
  if(!missing(tags) && !(missing(new.desc) && missing(new.labl) && missing(new.tags))) {
    nama <- names(x)
    desc <- attr(x, "desc")
    labl <- attr(x, "labl")
    valu <- attr(x, "valu")
    for(i in seq_along(tags))
    if(!is.na(m <- match(tags[i], nama))) {
      if(!missing(new.desc)) desc[m] <- new.desc[i]
      if(!missing(new.labl)) labl[m] <- new.labl[i]
      if(!missing(new.tags)) {
        names(x)[m] <- new.tags[i]
        if(tags[i] == valu)
          attr(x, "valu") <- new.tags[i]
      }
    }
    attr(x, "desc") <- desc
    attr(x, "labl") <- labl
  }
  if(!missing(new.dotnames))
    fvnames(x, ".") <- new.dotnames
  if(!missing(new.preferred)) {
    stopifnot(new.preferred %in% names(x))
    attr(x, "valu") <- new.preferred
  }
  if(!missing(new.formula))
    formula(x) <- new.formula
  return(x)
}

## common invocations to label a function like Kdot or Kcross
rebadge.as.crossfun <- function(x, main, sub=NULL, i, j) {
  i <- make.parseable(i)
  j <- make.parseable(j)
  if(is.null(sub)) {
    ylab <- substitute(main[i, j](r),
                       list(main=main, i=i, j=j))
    fname <- c(main, paste0("list", paren(paste(i, j, sep=","))))
    yexp <- substitute(main[list(i, j)](r),
                       list(main=main, i=i, j=j))
  } else {
    ylab <- substitute(main[sub, i, j](r),
                       list(main=main, sub=sub, i=i, j=j))
    fname <- c(main, paste0("list", paren(paste(sub, i, j, sep=","))))
    yexp <- substitute(main[list(sub, i, j)](r),
                       list(main=main, sub=sub, i=i, j=j))
  }
  y <- rebadge.fv(x, new.ylab=ylab, new.fname=fname, new.yexp=yexp)
  return(y)
}

rebadge.as.dotfun <- function(x, main, sub=NULL, i) {
  i <- make.parseable(i)
  if(is.null(sub)) {
    ylab <- substitute(main[i ~ dot](r),
                       list(main=main, i=i))
    fname <- c(main, paste0(i, "~symbol(\"\\267\")"))
    yexp <- substitute(main[i ~ symbol("\267")](r),
                       list(main=main, i=i))
  } else {
    ylab <- substitute(main[sub, i ~ dot](r),
                       list(main=main, sub=sub, i=i))
    fname <- c(main, paste0("list",
                            paren(paste0(sub, ",",
                                         i, "~symbol(\"\\267\")"))))
    yexp <- substitute(main[list(sub, i ~ symbol("\267"))](r),
                       list(main=main, sub=sub, i=i))
  }
  y <- rebadge.fv(x, new.ylab=ylab, new.fname=fname, new.yexp=yexp)
  return(y)
}

## even simpler wrapper for rebadge.fv
rename.fv <- function(x, fname, ylab, yexp=ylab) {
  stopifnot(is.fv(x))
  stopifnot(is.character(fname) && (length(fname) %in% 1:2))
  argu <- fvnames(x, ".x")
  if(missing(ylab) || is.null(ylab))
    ylab <- switch(length(fname),
                   substitute(fn(argu), list(fn=as.name(fname),
                                             argu=as.name(argu))),
                   substitute(fn[fsub](argu), list(fn=as.name(fname[1]),
                                                   fsub=as.name(fname[2]),
                                                   argu=as.name(argu))))
  if(missing(yexp) || is.null(yexp))
    yexp <- ylab
  y <- rebadge.fv(x, new.fname=fname, new.ylab=ylab, new.yexp=yexp)
  return(y)
}

## subset extraction operator
"[.fv" <-
  function(x, i, j, ..., drop=FALSE)
{
  igiven <- !missing(i)
  jgiven <- !missing(j)
  y <- as.data.frame(x)
  if(igiven && jgiven)
    z <- y[i, j, drop=drop]
  else if(igiven)
    z <- y[i, , drop=drop]
  else if(jgiven)
    z <- y[ , j, drop=drop]
  else z <- y

  ## return only the selected values as a data frame or vector.
  if(drop) return(z)

  if(!jgiven) 
    selected <- seq_len(ncol(x))
  else {
    nameindices <- seq_along(names(x))
    names(nameindices) <- names(x)
    selected <- as.vector(nameindices[j])
  }

  # validate choice of selected/dropped columns
  nama <- names(z)
  argu <- attr(x, "argu")
  if(!(argu %in% nama))
    stop(paste("The function argument", sQuote(argu), "must not be removed"))
  valu <- attr(x, "valu")
  if(!(valu %in% nama))
    stop(paste("The default column of function values",
               sQuote(valu), "must not be removed"))

  # if the plot formula involves explicit mention of dropped columns,
  # replace it by a generic formula
  fmla <- as.formula(attr(x, "fmla"))
  if(!all(variablesinformula(fmla) %in% nama)) 
    fmla <- as.formula(. ~ .x, env=environment(fmla))
  
  ## If range of argument was implicitly changed, adjust "alim"
  alim <- attr(x, "alim")
  rang <- range(z[[argu]])
  alim <- intersect.ranges(alim, rang, fatal=FALSE)

  result <- fv(z, argu=attr(x, "argu"),
               ylab=attr(x, "ylab"),
               valu=attr(x, "valu"),
               fmla=fmla,
               alim=alim,
               labl=attr(x, "labl")[selected],
               desc=attr(x, "desc")[selected],
               unitname=attr(x, "units"),
               fname=attr(x,"fname"),
               yexp=attr(x, "yexp"))
  
  ## carry over preferred names, if possible
  dotn <- fvnames(x, ".")
  fvnames(result, ".") <- dotn[dotn %in% colnames(result)]
  shad <- fvnames(x, ".s")
  if(length(shad) && all(shad %in% colnames(result)))
    fvnames(result, ".s") <- shad
  return(result)
}  

## Subset and column replacement methods
## to guard against deletion of columns

"[<-.fv" <- function(x, i, j, value) {
  if(!missing(j)) {
    ## check for alterations to structure of object
    if((is.character(j) && !all(j %in% colnames(x))) ||
       (is.numeric(j) && any(j > ncol(x))))
      stop("Use bind.fv to add new columns to an object of class fv")
    if(is.null(value) && missing(i)) {
      ## column(s) will be removed
      co <- seq_len(ncol(x))
      names(co) <- colnames(x)
      keepcol <- setdiff(co, co[j])
      return(x[ , keepcol, drop=FALSE])
    }
  }
  NextMethod("[<-")
}

"$<-.fv" <- function(x, name, value) {
  j <- which(colnames(x) == name)
  if(is.null(value)) {
    ## column will be removed
    if(length(j) != 0)
      return(x[, -j, drop=FALSE])
    return(x)
  }
  if(length(j) == 0) {
    ## new column
    df <- data.frame(1:nrow(x), value)[,-1L,drop=FALSE]
    colnames(df) <- name
    y <- bind.fv(x, df, desc=paste("Additional variable", sQuote(name)))
    return(y)
  }
  NextMethod("$<-")
}

## method for 'formula'

formula.fv <- function(x, ...) {
  attr(x, "fmla")
}

# new generic

"formula<-" <- function(x, ..., value) {
  UseMethod("formula<-")
}

"formula<-.fv" <- function(x, ..., value) {
  if(is.null(value))
    value <- paste(fvnames(x, ".y"), "~", fvnames(x, ".x"))
  else if(inherits(value, "formula")) {
    ## convert formula to string
    value <- flat.deparse(value)
  } else if(!is.character(value))
    stop("Assignment value should be a formula or a string")
  attr(x, "fmla") <- value
  return(x)
}

##   method for with()

  
with.fv <- function(data, expr, ..., fun=NULL, enclos=NULL) {
  if(any(names(list(...)) == "drop"))
    stop("Outdated argument 'drop' used in with.fv")
  cl <- short.deparse(sys.call())
  verifyclass(data, "fv")
  if(is.null(enclos)) 
    enclos <- parent.frame()
   ## convert syntactic expression to 'expression' object
#  e <- as.expression(substitute(expr))
  ## convert syntactic expression to call
  elang <- substitute(expr)
  ## map "." etc to names of columns of data
  datanames <- names(data)
  xname <- fvnames(data, ".x")
  yname <- fvnames(data, ".y")
  ux <- as.name(xname)
  uy <- as.name(yname)
  dnames <- intersect(datanames, fvnames(data, "."))
  ud <- as.call(lapply(c("cbind", dnames), as.name))
  anames <- intersect(datanames, fvnames(data, ".a"))
  ua <- as.call(lapply(c("cbind", anames), as.name))
  if(length(snames <- fvnames(data, ".s"))) {
    snames <- intersect(datanames, snames)
    us <- as.call(lapply(c("cbind", snames), as.name))
  } else us <- NULL
  expandelang <- eval(substitute(substitute(ee,
                                      list(.=ud, .x=ux, .y=uy, .s=us, .a=ua)),
                           list(ee=elang)))
  dont.complain.about(ua, ud, us, ux, uy)
  evars <- all.vars(expandelang)
  used.dotnames <- evars[evars %in% dnames]
  ## evaluate expression
  datadf <- as.data.frame(data)
  results <- eval(expandelang, as.list(datadf), enclos=enclos)
  ## --------------------
  ## commanded to return numerical values only?
  if(!is.null(fun) && !fun)
    return(results)

  if(!is.matrix(results) && !is.data.frame(results)) {
    ## result is a vector
    if(is.null(fun)) fun <- FALSE
    if(!fun || length(results) != nrow(datadf))
      return(results)
    results <- matrix(results, ncol=1)
  } else {
    ## result is a matrix or data frame
    if(is.null(fun)) fun <- TRUE
    if(!fun || nrow(results) != nrow(datadf))
      return(results)
  }
  ## result is a matrix or data frame of the right dimensions
  ## make a new fv object
  ## ensure columns of results have names
  if(is.null(colnames(results)))
    colnames(results) <- paste("col", seq_len(ncol(results)), sep="")
  resultnames <- colnames(results)
  ## get values of function argument
  xvalues <- datadf[[xname]]
  ## tack onto result matrix
  results <- cbind(xvalues, results)
  colnames(results) <- c(xname, resultnames)
  results <- data.frame(results)
  ## check for alteration of column names
  oldnames <- resultnames
  resultnames <- colnames(results)[-1L]
  if(any(resultnames != oldnames))
    warning("some column names were illegal and have been changed")
  ## determine mapping (if any) from columns of output to columns of input
  namemap <- match(colnames(results), names(datadf))
  okmap <- !is.na(namemap)
  ## Build up fv object
  ## decide which of the columns should be the preferred value
  newyname <- if(yname %in% resultnames) yname else resultnames[1L]
  ## construct default plot formula
  fmla <- flat.deparse(as.formula(paste(". ~", xname)))
  dotnames <- resultnames
  ## construct description strings
  desc <- character(ncol(results))
  desc[okmap] <- attr(data, "desc")[namemap[okmap]]
  desc[!okmap] <- paste("Computed value", resultnames[!okmap])
  ## function name (fname) and mathematical expression for function (yexp)
  oldyexp <- attr(data, "yexp")
  oldfname <- attr(data, "fname")
  if(is.null(oldyexp)) {
    fname <- cl
    yexp <- substitute(f(xname), list(f=as.name(fname), xname=as.name(xname)))
  } else {
    ## map 'cbind(....)' to "." for name of function only
    cb <- paste("cbind(",
                paste(used.dotnames, collapse=","),
                ")", sep="")
    compresselang <- gsub(cb, ".", flat.deparse(expandelang), fixed=TRUE)
    compresselang <- as.formula(paste(compresselang, "~1"))[[2L]]
    ## construct mapping using original function name
    labmap <- fvlabelmap(data, dot=TRUE)
    labmap[["."]] <- oldyexp
    yexp <- eval(substitute(substitute(ee, ff), 
                            list(ee=compresselang, ff=labmap)))
    labmap2 <- labmap
    labmap2[["."]] <- as.name(oldfname)
    fname <- eval(substitute(substitute(ee, ff), 
                             list(ee=compresselang,
                                  ff=labmap2)))
    fname <- paren(flat.deparse(fname))
  }
  ## construct mathematical labels
  mathlabl <- as.character(fvlegend(data, expandelang))
  mathlabl <- gsub("[[:space:]]+", " ", mathlabl)
  labl <- colnames(results)
  mathmap <- match(labl, used.dotnames)
  okmath <- !is.na(mathmap)
  labl[okmath] <- mathlabl[mathmap[okmath]]
  ## form fv object and return
  out <- fv(results, argu=xname, valu=newyname, labl=labl,
            desc=desc, alim=attr(data, "alim"), fmla=fmla,
            unitname=unitname(data), fname=fname, yexp=yexp, ylab=yexp)
  fvnames(out, ".") <- dotnames
  return(out)
}

## method for 'range'

range.fv <- local({

  getValues <- function(x) {
    xdat <- as.matrix(as.data.frame(x))
    yall <- fvnames(x, ".")
    vals <- xdat[, yall]
    return(as.vector(vals))
  }
  
  range.fv <- function(..., na.rm=TRUE, finite=na.rm) {
    aarg <- list(...)
    isfun <- sapply(aarg, is.fv)
    if(any(isfun)) 
      aarg[isfun] <- lapply(aarg[isfun], getValues)
    z <- do.call(range, append(aarg, list(na.rm=na.rm, finite=finite)))
    return(z)
  }

  range.fv
})

min.fv <- function(..., na.rm=TRUE, finite=na.rm) {
  range(..., na.rm=TRUE, finite=na.rm)[1L]
}

max.fv <- function(..., na.rm=TRUE, finite=na.rm) {
  range(..., na.rm=TRUE, finite=na.rm)[2L]
}

  
## stieltjes integration for fv objects

stieltjes <- function(f, M, ...) {
  ## stieltjes integral of f(x) dM(x)
  stopifnot(is.function(f))
  if(is.stepfun(M)) {
    envM <- environment(M)
    #' jump locations
    x <- get("x", envir=envM)
    #' values of integrand
    fx <- f(x, ...)
    #' jump amounts
    xx <- c(-Inf, (x[-1L] + x[-length(x)])/2, Inf)
    dM <- diff(M(xx))
    #' integrate f(x) dM(x)
    f.dM <- fx * dM
    result <- sum(f.dM[is.finite(f.dM)])
    return(list(result))
  } else if(is.fv(M)) {
    ## integration variable
    argu <- attr(M, "argu")
    x <- M[[argu]]
    ## values of integrand
    fx <- f(x, ...)
    ## estimates of measure
    valuenames <- names(M) [names(M) != argu]
    Mother <- as.data.frame(M)[, valuenames]
    Mother <- as.matrix(Mother, nrow=nrow(M))
    ## increments of measure
    dM <- apply(Mother, 2, diff)
    dM <- rbind(dM, 0)
    ## integrate f(x) dM(x)
    f.dM <- fx * dM
    f.dM[!is.finite(f.dM)] <- 0
    results <- colSums(f.dM)
    results <- as.list(results)
    names(results) <- valuenames
    return(results)
  } else stop("M must be an object of class fv or stepfun")
}

prefixfv <- function(x, tagprefix="", descprefix="", lablprefix=tagprefix,
                     whichtags=fvnames(x, "*")) {
  ## attach a prefix to fv information 
  stopifnot(is.fv(x))
  att <- attributes(x)
  relevant <- names(x) %in% whichtags
  oldtags <- names(x)[relevant]
  newtags <- paste(tagprefix, oldtags, sep="")
  newlabl <- paste(lablprefix, att$labl[relevant], sep="")
  newdesc <- paste(descprefix, att$desc[relevant])
  y <- rebadge.fv(x, tags=oldtags,
                  new.desc=newdesc,
                  new.labl=newlabl,
                  new.tags=newtags)
  return(y)
}

reconcile.fv <- local({

  reconcile.fv <- function(...) {
    ## reconcile several fv objects by finding the columns they share in common
    z <- list(...)
    if(!all(unlist(lapply(z, is.fv)))) {
      if(length(z) == 1 &&
         is.list(z[[1L]]) &&
         all(unlist(lapply(z[[1L]], is.fv))))
        z <- z[[1L]]
      else    
        stop("all arguments should be fv objects")
    }
    n <- length(z)
    if(n <= 1) return(z)
    ## find columns that are common to all estimates
    keepcolumns <- names(z[[1L]])
    keepvalues <- fvnames(z[[1L]], "*")
    for(i in 2:n) {
      keepcolumns <- intersect(keepcolumns, names(z[[i]]))
      keepvalues <- intersect(keepvalues, fvnames(z[[i]], "*"))
    }
    if(length(keepvalues) == 0)
      stop("cannot reconcile fv objects: they have no columns in common")
    ## determine name of the 'preferred' column
    prefs <- unlist(lapply(z, fvnames, a=".y"))
    prefskeep <- prefs[prefs %in% keepvalues]
    if(length(prefskeep) > 0) {
      ## pick the most popular
      chosen <- unique(prefskeep)[which.max(table(prefskeep))]
    } else {
      ## drat - pick a value arbitrarily
      chosen <- keepvalues[1L]
    }
    z <- lapply(z, rebadge.fv, new.preferred=chosen)
    z <- lapply(z, "[.fv", j=keepcolumns)
    ## also clip to the same r values
    rmax <- min(sapply(z, maxrval))
    z <- lapply(z, cliprmax, rmax=rmax)
    return(z)
  }

  maxrval <- function(x) { max(with(x, .x)) }
  cliprmax <- function(x, rmax) { x[ with(x, .x) <= rmax, ] }
  
  reconcile.fv
})

as.function.fv <- function(x, ..., value=".y", extrapolate=FALSE) {
  trap.extra.arguments(...)
  value.orig <- value
  ## extract function argument
  xx <- with(x, .x)
  ## extract all function values 
  yy <- as.data.frame(x)[, fvnames(x, "*"), drop=FALSE]
  ## determine which value(s) to supply
  if(!is.character(value))
    stop("value should be a string or vector specifying columns of x")
  if(!all(value %in% colnames(yy))) {
    expandvalue <- try(fvnames(x, value))
    if(!inherits(expandvalue, "try-error")) {
      value <- expandvalue
    } else stop("Unable to determine columns of x")
  }
  yy <- yy[,value, drop=FALSE]
  argname <- fvnames(x, ".x")
  ## determine extrapolation rule (1=NA, 2=most extreme value)
  stopifnot(is.logical(extrapolate))
  stopifnot(length(extrapolate) %in% 1:2)
  endrule <- 1 + extrapolate
  ## make function(s)
  if(length(value) == 1 && !identical(value.orig, "*")) {
    ## make a single 'approxfun' and return it
    f <- approxfun(xx, yy[,,drop=TRUE], rule=endrule)
    ## magic
    names(formals(f))[1L] <- argname
    body(f)[[4L]] <- as.name(argname)
  } else {
    ## make a list of 'approxfuns' with different function values
    funs <- lapply(yy, approxfun, x = xx, rule = endrule)
    ## return a function which selects the appropriate 'approxfun' and executes
    f <- function(xxxx, what=value) {
      what <- match.arg(what)
      funs[[what]](xxxx)
    }
    ## recast function definition
    ## ('any sufficiently advanced technology is
    ##   indistinguishable from magic' -- Arthur C. Clarke)
    formals(f)[[2L]] <- value
    names(formals(f))[1L] <- argname
    ##    body(f)[[3L]][[2L]] <- as.name(argname)
    body(f) <- eval(substitute(substitute(z,
                                          list(xxxx=as.name(argname))),
                               list(z=body(f))))
  }
  class(f) <- c("fvfun", class(f))
  attr(f, "fname") <- attr(x, "fname")
  attr(f, "yexp") <- attr(x, "yexp")
  return(f)
}

print.fvfun <- function(x, ...) {
  y <- args(x)
  yexp <- as.expression(attr(x, "yexp"))
  body(y) <- as.name(paste("Returns interpolated value of", yexp))
  print(y, ...)
  return(invisible(NULL))
}

findcbind <- function(root, depth=0, maxdepth=1000) {
  ## recursive search through a parse tree to find calls to 'cbind'
  if(depth > maxdepth) stop("Reached maximum depth")
  if(length(root) == 1) return(NULL)
  if(identical(as.name(root[[1L]]), as.name("cbind"))) return(list(numeric(0)))
  out <- NULL
  for(i in 2:length(root)) {
    di <- findcbind(root[[i]], depth+1, maxdepth)
    if(!is.null(di))
      out <- append(out, lapply(di, append, values=i, after=FALSE))
  }
  return(out)
}

.MathOpNames <- c("+", "-", "*", "/",
                  "^", "%%", "%/%",
                  "&", "|", "!",
                  "==", "!=", "<", "<=", ">=", ">")

distributecbind <- local({

  distributecbind <- function(x) {
    ## x is an expression involving a call to 'cbind'
    ## return a vector of expressions, each obtained by replacing 'cbind(...)'
    ## by one of its arguments in turn.
    stopifnot(typeof(x) == "expression")
    xlang <- x[[1L]]
    locations <- findcbind(xlang)
    if(length(locations) == 0)
      return(x)
    ## cbind might occur more than once
    ## check that the number of arguments is the same each time
    narg <- unique(sapply(locations, nargs.in.expr, e=xlang))
    if(length(narg) > 1) 
      return(NULL)
    out <- NULL
    if(narg > 0) {
      for(i in 1:narg) {
        ## make a version of the expression
        ## in which cbind() is replaced by its i'th argument
        fakexlang <- xlang
        for(loc in locations) {
          if(length(loc) > 0) {
            ## usual case: 'loc' is integer vector representing nested index
            cbindcall <- xlang[[loc]]
            ## extract i-th argument
            argi <- cbindcall[[i+1]]
            ## if argument is an expression, enclose it in parentheses
            if(length(argi) > 1 && paste(argi[[1L]]) %in% .MathOpNames)
              argi <- substitute((x), list(x=argi))
            ## replace cbind call by its i-th argument
            fakexlang[[loc]] <- argi
          } else {
            ## special case: 'loc' = integer(0) representing xlang itself
            cbindcall <- xlang
            ## extract i-th argument
            argi <- cbindcall[[i+1L]]
            ## replace cbind call by its i-th argument
            fakexlang <- cbindcall[[i+1L]]
          }
        }
        ## add to final expression
        out <- c(out, as.expression(fakexlang))
      }
    }
    return(out)
  }

  nargs.in.expr <- function(loc, e) {
    n <- if(length(loc) > 0) length(e[[loc]]) else length(e)
    return(n - 1L)
  }

  distributecbind
})

## Form a new 'fv' object as a ratio

ratfv <- function(df, numer, denom, ..., ratio=TRUE) {
  ## Determine y
  if(!missing(df) && !is.null(df)) {
    y <- fv(df, ...)
    num <- NULL
  } else {
    ## Compute numer/denom
    ## Numerator must be a data frame
    num <- fv(numer, ...)    
    ## Denominator may be a data frame or a constant
    force(denom)
    y <- eval.fv(num/denom)
    ## relabel
    y <- fv(as.data.frame(y), ...)
  }
  if(!ratio)
    return(y)
  if(is.null(num)) {
    ## Compute num = y * denom
    ## Denominator may be a data frame or a constant
    force(denom)
    num <- eval.fv(y * denom)
    ## ditch labels
    num <- fv(as.data.frame(num), ...)
  }
  ## make denominator an fv object
  if(is.data.frame(denom)) {
    den <- fv(denom, ...)
  } else {
    ## scalar
    check.1.real(denom, "Unless it is a data frame,")
    ## replicate it in all the data columns
    dendf <- as.data.frame(num)
    valuecols <- (names(num) != fvnames(num, ".x"))
    dendf[, valuecols] <- denom
    den <- fv(dendf, ...)
  } 
  ## tweak the descriptions
  ok <- (names(y) != fvnames(y, ".x"))
  attr(num, "desc")[ok] <- paste("numerator of",   attr(num, "desc")[ok])
  attr(den, "desc")[ok] <- paste("denominator of", attr(den, "desc")[ok])
  ## form ratio object
  y <- rat(y, num, den, check=FALSE)
  return(y)
}

## Tack new column(s) onto a ratio fv object

bind.ratfv <- function(x, numerator=NULL, denominator=NULL, 
                       labl = NULL, desc = NULL, preferred = NULL,
                       ratio=TRUE,
		       quotient=NULL) {
  if(ratio && !inherits(x, "rat"))
    stop("ratio=TRUE is set, but x has no ratio information", call.=FALSE)
  if(is.null(numerator) && !is.null(denominator) && !is.null(quotient))
    numerator <- quotient * denominator
  if(is.null(denominator) && inherits(numerator, "rat")) {
    ## extract numerator & denominator from ratio object
    both <- numerator
    denominator <- attr(both, "denominator")
    usenames <- fvnames(both, ".a")
    numerator   <- as.data.frame(both)[,usenames,drop=FALSE]
    denominator <- as.data.frame(denominator)[,usenames,drop=FALSE]
    ##  labels default to those of ratio object
    ma <- match(usenames, colnames(both))
    if(is.null(labl)) labl <- attr(both, "labl")[ma]
    if(is.null(desc)) desc <- attr(both, "desc")[ma]
  }
  # calculate ratio
  #    The argument 'quotient' is rarely needed 
  #    except to avoid 0/0 or to improve accuracy
  if(is.null(quotient))
    quotient <- numerator/denominator
    
  # bind new column to x   
  y <- bind.fv(x, quotient,
               labl=labl, desc=desc, preferred=preferred)
  if(!ratio)
    return(y)
    
  ## convert scalar denominator to data frame
  if(!is.data.frame(denominator)) {
    if(!is.numeric(denominator) || !is.vector(denominator))
      stop("Denominator should be a data frame or a numeric vector")
    nd <- length(denominator)
    if(nd != 1 && nd != nrow(x))
      stop("Denominator has wrong length")
    dvalue <- denominator
    denominator <- numerator
    denominator[] <- dvalue
  }
  ## Now fuse with x
  num <- attr(x, "numerator")
  den <- attr(x, "denominator")
  num <- bind.fv(num, numerator,
                 labl=labl, desc=paste("numerator of", desc),
                 preferred=preferred)
  den <- bind.fv(den, denominator,
                 labl=labl, desc=paste("denominator of", desc),
                 preferred=preferred)
  y <- rat(y, num, den, check=FALSE)
  return(y)
}

conform.ratfv <- function(x) {
  ## harmonise display properties in components of a ratio
  stopifnot(inherits(x, "rat"), is.fv(x))
  num <- attr(x, "numerator")
  den <- attr(x, "denominator")
  formula(num) <- formula(den) <- formula(x)
  fvnames(num, ".") <- fvnames(den, ".") <- fvnames(x, ".")
  unitname(num)     <- unitname(den)     <- unitname(x)
  attr(x, "numerator") <- num
  attr(x, "denominator") <- den
  return(x)
}

