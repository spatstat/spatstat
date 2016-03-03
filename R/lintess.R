#'
#'   lintess.R
#'
#'   Tessellations on a Linear Network
#'
#'   $Revision: 1.5 $   $Date: 2016/03/03 05:59:56 $
#'

lintess <- function(L, df) {
  verifyclass(L, "linnet")
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
      thedata <- paste("Data for segment", piece$seg[[1]])
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
  splat(nt, "tiles")
  return(invisible(NULL))
}

plot.lintess <- function(x, ..., main) {
  if(missing(main)) main <- short.deparse(substitute(x))
  plot(as.linfun(x), main=main, ...)
}

as.linfun.lintess <- function(X, ...) {
  L <- X$L
  df <- X$df
  f <- function(x, y, seg, tp) {
    result <- df$tile[integer(0)]
    for(i in seq_along(seg)) {
      tpi <- tp[i]
      segi <- seg[i]
      j <- which(df$seg == segi)
      kk <- which(df[j, "t0"] <= tpi & df[j, "t1"] >= tpi)
      result[i] <- if(length(kk) == 0) NA else df$tile[j[min(kk)]]
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
      jedge <- c(jedge, joined[-1])
    }
    nedge <- length(iedge)
    zz <- .C("cocoGraph",
             nv = as.integer(nbits),
             ne = as.integer(nedge), 
             ie = as.integer(iedge - 1L),
             je = as.integer(jedge - 1L),
             label = as.integer(integer(nbits)), 
             status = as.integer(integer(1)))
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
               data.frame(seg=c(seg[1], seg),
                          t0 = c(0, tp),
                          t1 = c(tp, 1),
                          from=NA_integer_,
                          to=NA_integer_))
    df$from[1] <- z$from[1]
    df$to[nrow(df)] <- z$to[1]
    return(df)
  }

  divide.linnet
})


