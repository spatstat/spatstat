#'
#'   lintess.R
#'
#'   Tessellations on a Linear Network
#'
#'   $Revision: 1.3 $   $Date: 2016/02/18 02:52:18 $
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

