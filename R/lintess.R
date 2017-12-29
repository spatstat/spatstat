#'
#'   lintess.R
#'
#'   Tessellations on a Linear Network
#'
#'   $Revision: 1.16 $   $Date: 2017/12/18 00:57:05 $
#'

lintess <- function(L, df) {
  verifyclass(L, "linnet")
  if(missing(df) || is.null(df)) {
    # tessellation consisting of a single tile
    ns <- nsegments(L)
    df <- data.frame(seg=seq_len(ns), t0=0, t1=1, tile=factor(1))
    out <- list(L=L, df=df)
    class(out) <- c("lintess", class(out))
    return(out)
  } 
  # validate 'df'
  stopifnot(is.data.frame(df))
  needed <- c("seg", "t0", "t1", "tile")
  if(any(bad <- is.na(match(needed, colnames(df)))))
    stop(paste(ngettext(sum(bad), "Column", "Columns"),
               commasep(sQuote(needed[bad])),
               "missing from data frame"),
         call.=FALSE)
  df$seg <- as.integer(df$seg)
  df$tile <- as.factor(df$tile)
  if(any(reversed <- with(df, t1 < t0)))
    df[reversed, c("t0", "t1")] <- df[reversed, c("t1", "t0")]
  with(df, {
    segU <- sort(unique(seg))
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
  out <- list(L=L, df=df)
  class(out) <- c("lintess", class(out))
  return(out)
}

print.lintess <- function(x, ...) {
  splat("Tessellation on a linear network")
  nt <- length(levels(x$df$tile))
  splat(nt, ngettext(nt, "tile", "tiles"))
  if(anyNA(x$df$tile)) splat("[An additional tile is labelled NA]")
  return(invisible(NULL))
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
  y <- list(nt=nt, nr=nr, lev=lev, seglen=seglen, tilelen=tilelen,
            hasna=hasna, nalen=nalen)
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
  })
  return(invisible(NULL))
}

plot.lintess <- function(x, ..., main, add=FALSE,
                         style=c("segments", "image"),
                         col=NULL) {
  if(missing(main)) main <- short.deparse(substitute(x))
  style <- match.arg(style)
  if(style == "image") {
    z <- plot(as.linfun(x), main=main, ..., add=add)
    return(invisible(z))
  }
  #' determine colour map
  df <- x$df
  lev <- levels(df$tile)
  if(is.null(col)) {
    col <- rainbow(length(lev))
    cmap <- colourmap(col, inputs=lev)
  } else if(inherits(col, "colourmap")) {
    cmap <- col
    col <- cmap(lev)
  } else if(is.colour(col)) {
    if(length(col) == 1) col <- rep(col, length(lev))
    if(length(col) != length(lev))
      stop(paste(length(col), "colours provided but",
                 length(lev), "colours needed"))
    cmap <- colourmap(col, inputs=lev)
  } else stop("col should be a vector of colours, or a colourmap object")
  #' determine segment coordinates
  L <- as.linnet(x)
  from <- L$from[df$seg]
  to   <- L$to[df$seg]
  V <- vertices(L)
  vx <- V$x
  vy <- V$y
  #' plot
  if(!add) plot(Frame(x), main=main, type="n")
  with(df,
       segments(
         vx[from] * (1-t0) + vx[to] * t0,
         vy[from] * (1-t0) + vy[to] * t0,
         vx[from] * (1-t1) + vx[to] * t1,
         vy[from] * (1-t1) + vy[to] * t1,
         col=col[as.integer(tile)],
         ...)
       )
  return(invisible(cmap))
}

as.owin.lintess <- function(W, ...) { as.owin(as.linnet(W), ...) }

Window.lintess <- function(X, ...) { as.owin(as.linnet(X)) }

domain.lintess <- as.linnet.lintess <- function(X, ...) { X$L }

as.linfun.lintess <- function(X, ..., values, navalue=NA) {
  L <- X$L
  df <- X$df
  if(missing(values) || is.null(values)) {
    rowvalues <- df$tile
  } else {
    if(length(values) != length(levels(df$tile)))
      stop("Length of 'values' should equal the number of tiles", call.=FALSE)
    rowvalues <- values[as.integer(df$tile)]    
  }
  f <- function(x, y, seg, tp) {
    result <- rowvalues[integer(0)]
    for(i in seq_along(seg)) {
      tpi <- tp[i]
      segi <- seg[i]
      j <- which(df$seg == segi)
      kk <- which(df[j, "t0"] <= tpi & df[j, "t1"] >= tpi)
      result[i] <- if(length(kk) == 0) navalue else rowvalues[j[min(kk)]]
    }
    return(result)
  }
  g <- linfun(f, L)
  return(g)
}

#'  Divide a linear network into tiles demarcated by
#' the points of a point pattern

divide.linnet <- local({
  
  divide.linnet <- function(X) {
    stopifnot(is.lpp(X))
    L <- as.linnet(X)
    coo <- coords(X)
    #' add identifiers of endpoints
    coo$from <- L$from[coo$seg]
    coo$to   <- L$to[coo$seg]
    #' group data by segment, sort by increasing 'tp'
    coo <- coo[with(coo, order(seg, tp)), , drop=FALSE]
    bits <- split(coo, coo$seg)
    #' expand as a sequence of intervals
    bits <- lapply(bits, expanddata)
    #' reassemble as data frame
    df <- Reduce(rbind, bits)
    #' find all undivided segments
    other <- setdiff(seq_len(nsegments(L)), unique(coo$seg))
    #' add a single line for each undivided segment
    if(length(other) > 0)
      df <- rbind(df, data.frame(seg=other, t0=0, t1=1,
                                 from=L$from[other], to=L$to[other]))
    #' We now have a tessellation 
    #' Sort again
    df <- df[with(df, order(seg, t0)), , drop=FALSE]
    #' Now identify connected components
    #' Two intervals are connected if they share an endpoint
    #' that is a vertex of the network.
    nvert <- nvertices(L)
    nbits <- nrow(df)
    iedge <- jedge <- integer(0)
    for(iv in seq_len(nvert)) {
      joined <- with(df, which(from == iv | to == iv))
      njoin <- length(joined)
      if(njoin > 1)
        iedge <- c(iedge, joined[-njoin])
      jedge <- c(jedge, joined[-1L])
    }
    nedge <- length(iedge)
    zz <- .C("cocoGraph",
             nv = as.integer(nbits),
             ne = as.integer(nedge), 
             ie = as.integer(iedge - 1L),
             je = as.integer(jedge - 1L),
             label = as.integer(integer(nbits)), 
             status = as.integer(integer(1L)),
             PACKAGE = "spatstat")
    if (zz$status != 0) 
      stop("Internal error: connectedness algorithm did not converge")
    lab <- zz$label + 1L
    lab <- as.integer(factor(lab))
    df <- df[,c("seg", "t0", "t1")]
    df$tile <- lab
    return(lintess(L, df))
  }

  expanddata <- function(z) {
    df <- with(z,
               data.frame(seg=c(seg[1L], seg),
                          t0 = c(0, tp),
                          t1 = c(tp, 1),
                          from=NA_integer_,
                          to=NA_integer_))
    df$from[1L] <- z$from[1L]
    df$to[nrow(df)] <- z$to[1L]
    return(df)
  }

  divide.linnet
})

