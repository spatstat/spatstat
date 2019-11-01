#'
#'   lintess.R
#'
#'   Tessellations on a Linear Network
#'
#'   $Revision: 1.41 $   $Date: 2019/09/17 07:23:27 $
#'

lintess <- function(L, df, marks=NULL) {
  verifyclass(L, "linnet")

  if(missing(df) || is.null(df)) {
    # tessellation consisting of a single tile
    ns <- nsegments(L)
    df <- data.frame(seg=seq_len(ns), t0=0, t1=1, tile=factor(1))
    return(lintess(L, df, marks))
  }
  
  #' validate 'df'
  stopifnot(is.data.frame(df))
  dfnames <- colnames(df)
  needed <- c("seg", "t0", "t1", "tile")
  if(any(bad <- is.na(match(needed, dfnames))))
    stop(paste(ngettext(sum(bad), "Column", "Columns"),
               commasep(sQuote(needed[bad])),
               "missing from data frame"),
         call.=FALSE)
  #' straighten out
  df <- df[, needed]
  df$seg <- as.integer(df$seg)
  df$tile <- as.factor(df$tile)
  if(any(reversed <- with(df, t1 < t0)))
    df[reversed, c("t0", "t1")] <- df[reversed, c("t1", "t0")]
  with(df, {
    segU <- sortunique(seg)
    segN <- seq_len(nsegments(L))
    if(length(omitted <- setdiff(segN, segU)) > 0)
      stop(paste(ngettext(length(omitted), "Segment", "Segments"),
                 commasep(omitted),
                 "omitted from data"),
           call.=FALSE)
    if(length(unknown <- setdiff(segU, segN)) > 0)
      stop(paste(ngettext(length(unknown), "Segment", "Segments"),
                 commasep(unknown),
                 ngettext(length(unknown), "do not", "does not"),
                 "exist in the network"),
           call.=FALSE)
    pieces <- split(df, seg)
    for(piece in pieces) {
      t0 <- piece$t0
      t1 <- piece$t1
      thedata <- paste("Data for segment", piece$seg[[1L]])
      if(!any(t0 == 0))
        stop(paste(thedata, "do not contain an entry with t0 = 0"),
             call.=FALSE)
      if(!any(t1 == 1))
        stop(paste(thedata, "do not contain an entry with t1 = 1"),
             call.=FALSE)
      if(any(t1 < 1 & is.na(match(t1, t0))) |
         any(t0 > 0 & is.na(match(t0, t1))))
        stop(paste(thedata, "are inconsistent"),
             call.=FALSE)
    }
  })
  #' validate marks
  if(!is.null(marks)) {
    marks <- as.data.frame(marks)
    nr <- nrow(marks)
    nt <- length(levels(df$tile))
    if(nr == 1L) {
      marks <- marks[rep(1, nt), , drop=FALSE]
      row.names(marks) <- 1:nt
      nr <- nt
    } else if(nr != nt) {
      stop(paste("Wrong number of",
                 ngettext(ncol(marks), "mark values:", "rows of mark values:"),
                 nr, "should be", nt),
           call.=FALSE)
    }
  }
  out <- list(L=L, df=df, marks=marks)
  class(out) <- c("lintess", class(out))
  return(out)
}

print.lintess <- function(x, ...) {
  splat("Tessellation on a linear network")
  nt <- length(levels(x$df$tile))
  splat(nt, ngettext(nt, "tile", "tiles"))
  if(anyNA(x$df$tile)) splat("[An additional tile is labelled NA]")
  if(!is.null(marx <- x$marks)) {
    mvt <- markvaluetype(marx)
    if(length(mvt) == 1) {
      splat("Tessellation has", mvt, "marks")
    } else {
      splat("Tessellation has mark variables",
            commasep(paste(colnames(marx), paren(mvt))))
    }
  }
  return(invisible(NULL))
}

nobjects.lintess <- function(x) {
  length(levels(x$df$tile))
}

tile.lengths <- function(x) {
  if(!inherits(x, "lintess"))
    stop("x should be a tessellation on a linear network (class 'lintess')",
         call.=FALSE)
  seglen <- lengths.psp(as.psp(x$L))
  df <- x$df
  df$fraglen <- with(df, seglen[seg] * (t1-t0))
  tilelen <- with(df, tapplysum(fraglen, list(tile)))
  return(tilelen)
}

tilenames.lintess <- function(x) {
 levels(x$df$tile) 
}

"tilenames<-.lintess" <- function(x, value) {
  levels(x$df$tile) <- value
  return(x)
}

marks.lintess <- function(x, ...) { x$marks }

"marks<-.lintess" <- function(x, ..., value) {
  if(!is.null(value)) {
    value <- as.data.frame(value)
    nt <- length(levels(x$df$tile))
    if(nrow(value) != nt)
      stop(paste("replacement value for marks has wrong length:",
                 nrow(value), "should be", nt),
           call.=FALSE)
    rownames(value) <- NULL
    if(ncol(value) == 1) colnames(value) <- "marks"
  }
  x$marks <- value
  return(x)
}

unmark.lintess <- function(X) {
  X$marks <- NULL
  return(X)
}

summary.lintess <- function(object, ...) {
  df <- object$df
  lev <- levels(df$tile)
  nt <- length(lev)
  nr <- nrow(df)
  seglen <- lengths.psp(as.psp(object$L))
  df$fraglen <- with(df, seglen[seg] * (t1-t0))
  tilelen <- with(df, tapplysum(fraglen, list(tile)))
  hasna <- anyNA(df$tile)
  nalen <- if(hasna) (sum(seglen) - sum(tilelen)) else 0
  marx <- object$marks
  if(!is.null(marx)) {
    mvt <- markvaluetype(marx)
    names(mvt) <- colnames(marx)
    marx <- summary(marx)
  } else mvt <- NULL
  y <- list(nt=nt, nr=nr, lev=lev, seglen=seglen, tilelen=tilelen,
            hasna=hasna, nalen=nalen, marx=marx, mvt=mvt)
  class(y) <- c("summary.lintess", class(y))
  return(y)
}

print.summary.lintess <- function(x, ...) {
  splat("Tessellation on a linear network")
  with(x, {
    splat(nt, "tiles")
    if(hasna) splat("[An additional tile is labelled NA]")
    if(nt <= 30) {
      splat("Tile labels:", paste(lev, collapse=" "))
      splat("Tile lengths:")
      print(signif(tilelen, 4))
    } else {
      splat("Tile lengths (summary):")
      print(summary(tilelen))
    }
    if(hasna) splat("Tile labelled NA has length", nalen)
    if(!is.null(marx)) {
      splat("Tessellation is marked")
      if(length(mvt) == 1) {
        splat("Marks are of type", sQuote(mvt))
      } else {
        splat("Mark variables:",
              commasep(paste(names(mvt), paren(unname(mvt)))))
      }
      splat("Summary of marks:")
      print(marx)
    }
  })
  return(invisible(NULL))
}

plot.lintess <- local({
  
  plot.lintess <- function(x, ..., main, add=FALSE,
                           style=c("colour", "width", "image"),
                           col=NULL,
                           values=marks(x),
                           ribbon=TRUE,
                           ribargs=list(),
                           multiplot=TRUE,
                           do.plot=TRUE
                           ) {
    if(missing(main)) main <- short.deparse(substitute(x))
    style <- match.arg(style)
    df <- x$df
    ntiles <- length(levels(df$tile))
    #' Associate 'values' with tiles
    if(markformat(values) == "hyperframe") 
      values <- as.data.frame(values) #' automatic warning
    switch(markformat(values),
           none = {
             #' no values assigned.
             #' default is tile name
             tn <- tilenames(x)
             values <- factor(tn, levels=tn)
           },
           vector = {
             #' vector of values.
             #' validate length of vector
             check.anyvector(values, ntiles, things="tiles")
           },
           dataframe = {
             #' data frame or matrix of values.
             values <- as.data.frame(values)
             if(nrow(values) != ntiles)
               stop(paste("Number of rows of values =", nrow(values),
                          "!=", ntiles, "= number of tiles"),
                    call.=FALSE)
               if(multiplot && ncol(values) > 1 && !add) {
                 #' Multiple Panel Plot
                 result <- multi.plot.lintess(x, ..., style=style, 
                                              main=main, do.plot=do.plot,
                                              ribbon=ribbon, ribargs=ribargs,
                                              col=col)
                 return(invisible(result))
               }
               if(ncol(values) > 1)
                 warning("Using only the first column of values")
                 values <- values[,1]
           },
           stop("Format of values is not understood")
           )
  
    #' Single Panel Plot
    if(style == "image") {
      z <- plot(as.linfun(x, values=values),
                main=main, ..., add=add, do.plot=do.plot,
                ribbon=ribbon, ribargs=ribargs, col=col)
      return(invisible(z))
    }
    #' convert to marked psp object
    L <- as.linnet(x)
    from <- L$from[df$seg]
    to   <- L$to[df$seg]
    V <- vertices(L)
    vx <- V$x
    vy <- V$y
    segdata <- with(df, list(x0=vx[from] * (1-t0) + vx[to] * t0,
                             y0=vy[from] * (1-t0) + vy[to] * t0,
                             x1=vx[from] * (1-t1) + vx[to] * t1,
                             y1=vy[from] * (1-t1) + vy[to] * t1,
                             marks=values[as.integer(tile)]))
    S <- as.psp(segdata, window=Window(L))
    cmap <- plot(S, style=style, add=add, do.plot=do.plot, main=main, 
                 ribbon=ribbon, ribargs=ribargs, col=col, ...)
    return(invisible(cmap))
  }

  multi.plot.lintess <- function(x, ...,
                                 zlim=NULL, col=NULL, equal.ribbon=FALSE) {
    if(equal.ribbon && is.null(zlim) && !inherits(col, "colourmap"))
      zlim <- range(marks(x))
    if(!is.null(zlim)) {
      result <- plot(unstack(x), ..., zlim=zlim, col=col)
    } else {
      result <- plot(unstack(x), ..., col=col)
    }
    return(invisible(result))
  }
  
  plot.lintess
})

  
unstack.lintess <- function(x, ...) {
  marx <- marks(x)
  if(is.null(marx) || is.null(dim(marx)) || ncol(marx) <= 1)
    return(solist(x))
  ux <- unmark(x)
  y <- solapply(as.list(marx), setmarks, x=ux)
  return(y)
} 

as.owin.lintess <- function(W, ...) { as.owin(as.linnet(W), ...) }

Window.lintess <- function(X, ...) { as.owin(as.linnet(X)) }

domain.lintess <- as.linnet.lintess <- function(X, ...) { X$L }

as.data.frame.lintess <- function(x, ...) {
  df <- x$df
  if(!is.null(marx <- marks(x))) {
    marx <- as.data.frame(marx)
    if(ncol(marx) == 1L) colnames(marx) <- "marks"
    marx <- marx[as.integer(df$tile), , drop=FALSE]
    df <- cbind(df, marx)
  }
  df <- as.data.frame(df, ...)
  return(df)
}

lineartileindex <- function(seg, tp, Z,
                            method=c("encode", "C", "interpreted")) {
  method <- match.arg(method)
  df <- if(inherits(Z, "lintess")) Z$df else
        if(is.data.frame(Z)) Z else stop("Format of Z is unrecognised")
  switch(method,
         interpreted = {
           n <- length(seg)
           #' extract tessellation data
           tilenames <- levels(df$tile)
           answer <- factor(rep(NA_integer_, n),
                            levels=seq_along(tilenames),
                            labels=tilenames)
           for(i in seq_along(seg)) {
             tpi <- tp[i]
             segi <- seg[i]
             j <- which(df$seg == segi)
             kk <- which(df[j, "t0"] <= tpi & df[j, "t1"] >= tpi)
             answer[i] <- if(length(kk) == 0) NA else df[j[min(kk)], "tile"]
           }
         },
         encode = {
           #' encode locations as numeric
           loc <- seg - 1 + tp
           #' extract tessellation data and sort them
           df <- df[order(df$seg, df$t0), , drop=FALSE]
           m <- nrow(df)
           #' encode breakpoints as numeric
           bks <- with(df, c(seg - 1 + t0, seg[m]))
           #' which piece contains each query point
           jj <- findInterval(loc, bks,
                              left.open=TRUE, all.inside=TRUE,
                              rightmost.closed=TRUE)
           answer <- df$tile[jj]
         },
         C = {
           #' sort query data
           oo <- order(seg, tp)
           seg <- seg[oo]
           tp  <- tp[oo]
           n <- length(seg)
           #' extract tessellation data and sort them
           df <- df[order(df$seg, df$t0), , drop=FALSE]
           m <- nrow(df)
           #' handle factor 
           dftile <- df$tile
           tilecode <- as.integer(dftile)
           tilenames <- levels(dftile)
           #' launch
           z <- .C("lintileindex",
                   n      = as.integer(n),
                   seg    = as.integer(seg),
                   tp     = as.double(tp),
                   dfn    = as.integer(m),
                   dfseg  = as.integer(df$seg),
                   dft0   = as.double(df$t0),
                   dft1   = as.double(df$t1),
                   dftile = as.integer(tilecode),
                   answer = as.integer(integer(n)),
                   PACKAGE="spatstat")
           z <- z$answer
           z[z == 0] <- NA
           answer <- integer(n)
           answer[oo] <- z
           answer <- factor(answer,
                            levels=seq_along(tilenames),
                            labels=tilenames)
           })
  return(answer)
}

as.linfun.lintess <- local({

  as.linfun.lintess <- function(X, ..., values=marks(X), navalue=NA) {
    Xdf <- X$df
    tilenames <- levels(Xdf$tile)
    value.is.tile <- is.null(values)
    if(value.is.tile) {
      tilevalues <- factor(tilenames, levels=tilenames)
    } else {
      if(!is.null(dim(values))) values <- values[,1]
      if(length(values) != length(tilenames))
        stop("Length of 'values' should equal the number of tiles", call.=FALSE)
      tilevalues <- values
    }
    f <- function(x, y, seg, tp) {
      k <- as.integer(lineartileindex(seg, tp, Xdf))
      if(!anyNA(k)) {
        result <- tilevalues[k]
      } else {
        ok <- !is.na(k)
        result <- rep(navalue, length(seg))
        result[ok] <- tilevalues[k[ok]]
      }
      return(result)
    }
    g <- linfun(f, X$L)
    attr(g, "explain") <- uitleggen
    return(g)
  }

  uitleggen <- function(x, ...) {
    envf <- environment(attr(x, "f"))
    Xdf <- get("Xdf", envir=envf)
    value.is.tile <- get("value.is.tile", envir=envf) %orifnull% FALSE
    if(value.is.tile) {
      valuename <- "the tile index"
    } else {
      tilevalues <- get("tilevalues", envir=envf)
      valuetype <- typeof(tilevalues)
      valuename <- paste("a value", paren(sQuote(valuetype)),
                         "associated with each tile")
    }
    splat("Function on a network, associated with a tessellation,")
    splat("\treturns", valuename)
    nt <- length(levels(Xdf$tile))
    splat("Tessellation has", nt, ngettext(nt, "tile", "tiles"))
    splat("Function domain:")
    print(as.linnet(x))
    return(invisible(NULL))
  }

  as.linfun.lintess
})

