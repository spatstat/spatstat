#
#  hyperframe.R
#
# $Revision: 1.52 $  $Date: 2014/10/24 00:22:30 $
#

hyperframe <- function(...,
                       row.names=NULL, check.rows=FALSE, check.names=TRUE,
                       stringsAsFactors=default.stringsAsFactors()) {
  aarg <- list(...)
  nama <- names(aarg)

  # number of columns (= variables)
  nvars <- length(aarg)
  
  if(nvars == 0) {
    # zero columns - return
    result <- list(nvars=0,
                   ncases=0,
                   vname=character(0),
                   vtype=factor(,
                     levels=c("dfcolumn","hypercolumn","hyperatom")),
                   vclass=character(0),
                   df=data.frame(),
                   hyperatoms=list(),
                   hypercolumns=list())
    class(result) <- c("hyperframe", class(result))
    return(result)
  }

  # check column names
  if(is.null(nama))
    nama <- paste("V", 1:nvars, sep="")
  else if(any(unnamed <- (nama == ""))) 
    nama[unnamed] <- paste("V", seq_len(sum(unnamed)), sep="")
  nama <- make.names(nama, unique=TRUE)
  names(aarg) <- nama
  
  # Each argument must be either
  #    - a vector suitable as a column in a data frame
  #    - a list of objects of the same class
  #    - a single object of some class
  
  is.dfcolumn <- function(x) {
    is.atomic(x) && (is.vector(x) || is.factor(x))
  }
  is.hypercolumn <- function(x) {
    if(!is.list(x))
      return(FALSE)
    if(length(x) <= 1)
      return(TRUE)
    cla <- class(x[[1]])
    all(sapply(x, function(xi,cla) { identical(class(xi), cla) }, cla=cla))
  }
  dfcolumns    <- sapply(aarg, is.dfcolumn)
  hypercolumns <- sapply(aarg, is.hypercolumn)
  hyperatoms   <- !(dfcolumns | hypercolumns)

  # Determine number of rows (= cases) 
  columns <- dfcolumns | hypercolumns
  if(!any(columns)) 
    ncases <- 1
  else {
    heights <- rep.int(1, nvars)
    heights[columns] <-  unlist(lapply(aarg[columns], length))
    u <- unique(heights)
    if(length(u) > 1) {
      u <- u[u != 1]
      if(length(u) > 1)
        stop(paste("Column lengths are inconsistent:",
                   paste(u, collapse=",")))
    }
    ncases <- u
    if(ncases > 1 && all(heights[dfcolumns] == 1)) {
      # force the data frame to have 'ncases' rows
      aarg[dfcolumns] <- lapply(aarg[dfcolumns], rep, ncases)
      heights[dfcolumns] <- ncases
    }
    if(any(stubs <- hypercolumns & (heights != ncases))) {
      ## hypercolumns of height 1 should be hyperatoms
      aarg[stubs] <- lapply(aarg[stubs], "[[", 1)
      hypercolumns[stubs] <- FALSE
      hyperatoms[stubs] <- TRUE
    }
  }

  # Collect the data frame columns into a data frame
  if(!any(dfcolumns))
    df <- as.data.frame(matrix(, ncases, 0), row.names=row.names)
  else {
    df <- do.call("data.frame", append(aarg[dfcolumns],
                                     list(row.names=row.names,
                                          check.rows=check.rows,
                                          check.names=check.names,
                                          stringsAsFactors=stringsAsFactors)))
    names(df) <- nama[dfcolumns]
  }

  # Storage type of each variable
  vtype <- character(nvars)
  vtype[dfcolumns] <- "dfcolumn"
  vtype[hypercolumns] <- "hypercolumn"
  vtype[hyperatoms] <- "hyperatom"
  vtype=factor(vtype, levels=c("dfcolumn","hypercolumn","hyperatom"))

  # Class of each variable
  class1 <- function(x) { class(x)[1] }
  vclass <- character(nvars)
  if(any(dfcolumns))
    vclass[dfcolumns] <- unlist(lapply(as.list(df), class1))
  if(any(hyperatoms))
    vclass[hyperatoms] <- unlist(lapply(aarg[hyperatoms], class1))
  if(any(hypercolumns))
    vclass[hypercolumns] <- unlist(lapply(aarg[hypercolumns],
                                          function(x) { class1(x[[1]]) }))
  
  # Put the result together
  result <- list(nvars=nvars,
                 ncases=ncases,
                 vname=nama,
                 vtype=vtype,
                 vclass=vclass,
                 df=df,
                 hyperatoms=aarg[hyperatoms],
                 hypercolumns=aarg[hypercolumns])
  
    class(result) <- c("hyperframe", class(result))
    return(result)
}

is.hyperframe <- function(x) inherits(x, "hyperframe")

print.hyperframe <- function(x, ...) {
  ux <- unclass(x)
  nvars <- ux$nvars
  ncases <- ux$ncases
  if(nvars * ncases == 0) {
    cat(paste("NULL hyperframe with", ncases,
              ngettext(ncases, "row (=case)", "rows (=cases)"),
              "and", nvars,
              ngettext(nvars, "column (=variable)", "columns (=variables)"),
              "\n"))
  } else {
    cat("Hyperframe:\n")
    print(as.data.frame(x, discard=FALSE), ...)
  }
  return(invisible(NULL))
}

dim.hyperframe <- function(x) {
  with(unclass(x), c(ncases, nvars))
}

summary.hyperframe <- function(object, ..., brief=FALSE) {
  x <- unclass(object)
  y <- list(
            nvars = x$nvars,
            ncases = x$ncases,
            dim = c(x$ncases, x$nvars),
            typeframe = data.frame(VariableName=x$vname, Class=x$vclass),
            storage = x$vtype,
            col.names = x$vname)
  classes <- x$vclass
  names(classes) <- x$vname
  y$classes <- classes
  # Ordinary data frame columns
  df <- x$df
  y$dfnames <- names(df)
  y$df <- if(length(df) > 0 && !brief) summary(df) else NULL
  y$row.names <- row.names(df)
  class(y) <- c("summary.hyperframe", class(y))
  return(y)
}

print.summary.hyperframe <- function(x, ...) {
  nvars <- x$nvars
  ncases <- x$ncases
  cat(paste(if(nvars * ncases == 0) "NULL" else NULL,
            "hyperframe with", ncases,
            ngettext(ncases, "row (=case)", "rows (=cases)"),
            "and", nvars,
            ngettext(nvars, "column (=variable)", "columns (=variables)"),
            "\n"))
  if(nvars == 0)
    return(invisible(NULL))
  # Variable names and types
  print(x$typeframe)
  # Ordinary data frame columns
  if(!is.null(x$df)) {
    cat("Summary of data frame columns:\n")
    print(x$df, ...)
  }
  return(invisible(NULL))
}

names.hyperframe <- function(x) { unclass(x)$vname }

"names<-.hyperframe" <- function(x, value) {
  x <- unclass(x)
  stopifnot(is.character(value))
  value <- make.names(value)
  if(length(value) != x$nvars)
    stop("Incorrect length for vector of names")
  vtype <- x$vtype
  names(x$df)           <- value[vtype == "dfcolumn"]
  names(x$hyperatoms)   <- value[vtype == "hyperatom"]
  names(x$hypercolumns) <- value[vtype == "hypercolumn"]
  x$vname <- value
  class(x) <- c("hyperframe", class(x))
  return(x)
}

row.names.hyperframe <- function(x) {
  return(row.names(unclass(x)$df))
}

"row.names<-.hyperframe" <- function(x, value) {
  y <- unclass(x)
  df <- y$df
  row.names(df) <- value
  y$df <- df
  class(y) <- c("hyperframe", class(y))
  return(y)
}


## conversion to hyperframe

as.hyperframe <- function(x, ...) {
  UseMethod("as.hyperframe")
}

as.hyperframe.hyperframe <- function(x, ...) {
  return(x)
}

as.hyperframe.data.frame <- function(x, ..., stringsAsFactors=FALSE) {
  xlist <- if(missing(x)) NULL else as.list(x)
  do.call("hyperframe",
          resolve.defaults(
                           xlist,
                           list(...),
                           list(row.names=rownames(x),
                                stringsAsFactors=stringsAsFactors),
                           .StripNull=TRUE))
}

as.hyperframe.anylist <- 
as.hyperframe.listof <- function(x, ...) {
  if(!missing(x)) {
    xname <- sensiblevarname(short.deparse(substitute(x)), "x")
    xlist <- list(x)
    names(xlist) <- xname
  } else xlist <- NULL
  do.call("hyperframe",
          resolve.defaults(
                           xlist,
                           list(...),
                           list(row.names=rownames(x)),
                           .StripNull=TRUE))
}

as.hyperframe.default <- function(x, ...) {
  as.hyperframe(as.data.frame(x, ...))
}

#### conversion to other types

as.data.frame.hyperframe <- function(x, row.names = NULL,
                                     optional = FALSE, ...,
                                     discard=TRUE, warn=TRUE) {
  ux <- unclass(x)
  if(is.null(row.names))
    row.names <- row.names(ux$df)
  vtype <- ux$vtype
  vclass <- ux$vclass
  dfcol <- (vtype == "dfcolumn")
  if(discard) { 
    nhyper <- sum(!dfcol)
    if(nhyper > 0 && warn)
      warning(paste(nhyper, 
                    ngettext(nhyper, "variable", "variables"),
                    "discarded in conversion to data frame"))
    df <- as.data.frame(ux$df, row.names=row.names, optional=optional, ...)
  } else {
    lx <- as.list(x)
    nrows <- ux$ncases
    vclassstring <- paren(vclass)
    if(any(!dfcol)) 
      lx[!dfcol] <- lapply(as.list(vclassstring[!dfcol]),
                           rep.int, times=nrows)
    df <- do.call("data.frame", append(lx, list(row.names=row.names)))
    colnames(df) <- ux$vname
  }
  return(df)
}

as.list.hyperframe <- function(x, ...) {
  ux <- unclass(x)
  nama <- ux$vname
  names(nama) <- nama
  out <- lapply(nama, function(nam, x) { x[, nam, drop=TRUE] }, x=x)
  if(ux$ncases == 1)
    out <- lapply(out, listof)
  out <- lapply(out, "names<-", value=row.names(x))
  return(out)
}

# evaluation

eval.hyper <- function(e, h, simplify=TRUE, ee=NULL) {
  .Deprecated("with.hyperframe", package="spatstat")
  if(is.null(ee))
    ee <- as.expression(substitute(e))
  with.hyperframe(h, simplify=simplify, ee=ee)
}

with.hyperframe <- function(data, expr, ..., simplify=TRUE, ee=NULL,
                            enclos=NULL) {
  if(!inherits(data, "hyperframe"))
    stop("data must be a hyperframe")
  if(is.null(ee))
    ee <- as.expression(substitute(expr))
  if(is.null(enclos))
    enclos <- parent.frame()
  n <- nrow(data)
  out <- vector(mode="list", length=n)
  datalist <- as.list(data)
  for(i in 1:n) {
    rowi <- lapply(datalist, "[[", i)  # ensures the result is always a list
    outi <- eval(ee, rowi, enclos)
    if(!is.null(outi))
      out[[i]] <- outi
  }
  names(out) <- row.names(data)
  if(simplify && all(unlist(lapply(out, is.vector)))) {
    # if all results are atomic vectors of equal length,
    # return a matrix or vector.
    lenfs <- unlist(lapply(out, length))
    if(all(unlist(lapply(out, is.atomic))) &&
            length(unique(lenfs)) == 1) {
      out <- t(as.matrix(as.data.frame(out)))
      row.names(out) <- row.names(data)
      out <- out[,,drop=TRUE]
      return(out)
    }
  }
  out <- hyperframe(result=out, row.names=row.names(data))$result
  return(out)
}

cbind.hyperframe <- function(...) {
  aarg <- list(...)
  narg <- length(aarg)
  if(narg == 0)
    return(hyperframe())
  namarg <- names(aarg)
  if(is.null(namarg))
    namarg <- rep.int("", narg)
  ishyper <- unlist(lapply(aarg, inherits, what="hyperframe"))
  isdf <- unlist(lapply(aarg, inherits, what="data.frame"))
  columns <- list()
  for(i in 1:narg) {
    if(ishyper[i] || isdf[i]){
      if(ncol(aarg[[i]]) > 0) {
        newcolumns <- as.list(aarg[[i]])
        if(namarg[i] != "")
          names(newcolumns) <- paste(namarg[i], ".", names(newcolumns), sep="")
        columns <- append(columns, newcolumns)
      }
    } else {
      nextcolumn <- list(aarg[[i]])
      names(nextcolumn) <- namarg[i]
      columns <- append(columns, nextcolumn)
    }
  }
  result <- do.call("hyperframe", columns)
  return(result)
}

rbind.hyperframe <- function(...) {
  argh <- list(...)
  if(length(argh) == 0)
    return(NULL)
  # convert them all to hyperframes
  argh <- lapply(argh, as.hyperframe)
  #
  nargh <- length(argh)
  if(nargh == 1)
    return(argh[[1]])
  # check for compatibility of dimensions & names
  dfs <- lapply(argh, as.data.frame, discard=FALSE)
  dfall <- do.call(rbind, dfs)
  # check that data frame columns also match
  dfs0 <- lapply(argh, as.data.frame, discard=TRUE, warn=FALSE)
  df0all <- do.call(rbind, dfs0)
  # assemble data
  rslt <- list()
  nam <- names(dfall) 
  nam0 <- names(df0all)
  for(k in seq_along(nam)) {
    nama <- nam[k]
    if(nama %in% nam0) {
      # data frame column: already made
      rslt[[k]] <- dfall[,k]
    } else {
      # hypercolumns or hyperatoms: extract them
      hdata <- lapply(argh,
                      function(x,nama) { x[, nama, drop=FALSE] },
                      nama=nama)
      hdata <- lapply(lapply(hdata, as.list), getElement, name=nama)
      # append them
      hh <- hdata[[1]]
      for(j in 2:nargh) {
        hh <- append(hh, hdata[[j]])
      }
      rslt[[k]] <- hh
    }
  }
  # make hyperframe
  names(rslt) <- nam
  out <- do.call(hyperframe, append(rslt, list(stringsAsFactors=FALSE)))
  return(out)
}

plot.hyperframe <-
  function(x, e, ..., main, arrange=TRUE,
           nrows=NULL, ncols=NULL,
           parargs=list(mar=c(1,1,3,1) * marsize),
           marsize=0.1) {
  xname <- short.deparse(substitute(x))
  main <- if(!missing(main)) main else xname

  if(missing(e)) {
    # default: plot first column that contains objects
    ok <- (summary(x)$storage %in% c("hypercolumn", "hyperatom"))
    if(any(ok)) {
      j <- min(which(ok))
      x <- x[,j, drop=TRUE]
      x <- as.listof(x)
      plot(x, ..., main=main, arrange=arrange, nrows=nrows, ncols=ncols)
      return(invisible(NULL))
    } else {
      # hyperframe does not contain any objects
      # invoke plot.data.frame
      x <- as.data.frame(x)
      plot(x, ..., main=main)
      return(invisible(NULL))
    }
  }

  if(!is.language(e))
    stop(paste("Argument e should be a call or an expression;",
               "use quote(...) or expression(...)"))
  ee <- as.expression(e)

  if(!arrange) {
    # No arrangement specified: just evaluate the plot expression 'nr' times
    with(x, ee=ee)
    return(invisible(NULL))
  }

  # Arrangement
  # Decide whether to plot a main header
  banner <- (sum(nchar(as.character(main))) > 0)
  if(length(main) > 1)
    main <- paste(main, collapse="\n")
  nlines <- if(!is.character(main)) 1 else length(unlist(strsplit(main, "\n")))
  # determine arrangement of plots
  # arrange like mfrow(nrows, ncols) plus a banner at the top
  n <- summary(x)$ncases
  if(is.null(nrows) && is.null(ncols)) {
    nrows <- as.integer(floor(sqrt(n)))
    ncols <- as.integer(ceiling(n/nrows))
  } else if(!is.null(nrows) && is.null(ncols))
    ncols <- as.integer(ceiling(n/nrows))
  else if(is.null(nrows) && !is.null(ncols))
    nrows <- as.integer(ceiling(n/ncols))
  else stopifnot(nrows * ncols >= length(x))
  nblank <- ncols * nrows - n
  # declare layout
  mat <- matrix(c(seq_len(n), numeric(nblank)),
                byrow=TRUE, ncol=ncols, nrow=nrows)
  heights <- rep.int(1, nrows)
  if(banner) {
    # Increment existing panel numbers
    # New panel 1 is the banner
    panels <- (mat > 0)
    mat[panels] <- mat[panels] + 1
    mat <- rbind(rep.int(1,ncols), mat)
    heights <- c(0.1 * (1 + nlines), heights)
  }
  # initialise plot
  layout(mat, heights=heights)
  # plot banner
  if(banner) {
    opa <- par(mar=rep.int(0,4), xpd=TRUE)
    plot(numeric(0),numeric(0),type="n",ann=FALSE,axes=FALSE,
         xlim=c(-1,1),ylim=c(-1,1))
    cex <- resolve.defaults(list(...), list(cex.title=2))$cex.title
    text(0,0,main, cex=cex)
  }
  # plot panels
  npa <- do.call("par", parargs)
  if(!banner) opa <- npa
  with(x, ee=ee)
  # revert
  layout(1)
  par(opa)
  return(invisible(NULL))
}


str.hyperframe <- function(object, ...) {
  d <- dim(object)
  x <- unclass(object)
  argh <- resolve.defaults(list(...), list(nest.lev=0, indent.str="  .."))
  cat(paste("'hyperframe':\t",
            d[1], ngettext(d[1], "row", "rows"),
            "and",
            d[2], ngettext(d[2], "column", "columns"),
            "\n"))
  nr <- d[1]
  nc <- d[2]
  if(nc > 0) {
    vname <- x$vname
    vclass <- x$vclass
    vtype  <- as.character(x$vtype)
    indentstring <- with(argh, paste(rep.int(indent.str, nest.lev), collapse=""))
    for(j in 1:nc) {
      tag <- paste("$", vname[j])
      switch(vtype[j],
             dfcolumn={
               desc <- vclass[j]
               if(nr > 0) {
                 vals <- object[1:min(nr,3),j,drop=TRUE]
                 vals <- paste(paste(format(vals), collapse=" "), "...")
               } else vals <- ""
             },
             hypercolumn=,
             hyperatom={
               desc <- "objects of class"
               vals <- vclass[j]
             })
      cat(paste(paste(indentstring, tag, sep=""),
                ":", desc, vals, "\n"))
    }
  }
  return(invisible(NULL))
}
