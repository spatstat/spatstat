#
#    summary.im.R
#
#    summary() method for class "im"
#
#    $Revision: 1.18 $   $Date: 2014/11/11 03:07:22 $
#
#    summary.im()
#    print.summary.im()
#    print.im()
#
summary.im <- function(object, ...) {
  verifyclass(object, "im")

  x <- object

  y <- unclass(x)[c("dim", "xstep", "ystep")]
  pixelarea <- y$xstep * y$ystep

  # extract image values
  v <- x$v
  inside <- !is.na(v)
  v <- v[inside]

  # type of values?
  y$type <- x$type
  
  # factor-valued?
  lev <- levels(x)
  if(!is.null(lev) && !is.factor(v))
    v <- factor(v, levels=seq_along(lev), labels=lev)

  switch(x$type,
         integer=,
         real={
           y$integral <- sum(v) * pixelarea
           y$mean <- mean(v)
           y$range <- range(v)
           y$min <- y$range[1]  
           y$max <- y$range[2]
         },
         factor={
           y$levels <- lev
           y$table <- table(v, dnn="")
         },
         complex={
           y$integral <- sum(v) * pixelarea
           y$mean <- mean(v)
           rr <- range(Re(v))
           y$Re <- list(range=rr, min=rr[1], max=rr[2])
           ri <- range(Im(v))
           y$Im <- list(range=ri, min=ri[1], max=ri[2])
         },
         {
           # another unknown type
           pixelvalues <- v
           y$summary <- summary(pixelvalues)
         })
    
  # summarise pixel raster
  win <- as.owin(x)
  y$window <- summary.owin(win)

  y$fullgrid <- (rescue.rectangle(win)$type == "rectangle")

  y$units <- unitname(x)
  
  class(y) <- "summary.im"
  return(y)
}

print.summary.im <- function(x, ...) {
  verifyclass(x, "summary.im")
  cat(paste0(x$type, "-valued"), "pixel image", fill=TRUE)
  unitinfo <- summary(x$units)
  pluralunits <- unitinfo$plural
  sigdig <- 5
  di <- x$dim
  win <- x$window
  cat(di[1], "x", di[2], "pixel array (ny, nx)", fill=TRUE)
  cat("enclosing rectangle:",
      prange(signif(x$window$xrange, sigdig)),
      "x",
      prange(signif(x$window$yrange, sigdig)),
      unitinfo$plural,
      unitinfo$explain,
      fill=TRUE)
  cat("dimensions of each pixel:",
      signif(x$xstep, 3), "x", signif(x$ystep, 3),
      pluralunits,
      fill=TRUE)
  if(!is.null(explain <- unitinfo$explain))
    cat(paste(explain, "\n"))
  if(x$fullgrid) {
    cat("Image is", "defined on the",  "full rectangular", "grid", fill=TRUE)
    whatpart <- "Frame"
  } else {
    cat("Image is", "defined on", "a subset of", "the rectangular", "grid",
        fill=TRUE)
    whatpart <- "Subset"
  }
  cat(whatpart, "area =", win$area, "square", pluralunits, fill=TRUE)
  if(x$fullgrid) cat("Pixel values:\n") else
                 cat("Pixel values", "(inside window):", fill=TRUE)
  switch(x$type,
         integer=,
         real={
           cat("\trange =", prange(signif(x$range, sigdig)), fill=TRUE)
           cat("\tintegral =", signif(x$integral, sigdig), fill=TRUE)
           cat("\tmean =", signif(x$mean, sigdig), fill=TRUE)
         },
         factor={
           print(x$table)
         },
         complex={
           cat("\trange: Real",
               prange(signif(x$Re$range, sigdig)),
               "Imaginary",
               prange(signif(x$Im$range, sigdig)),
               fill=TRUE)
           cat("\tintegral =", signif(x$integral, sigdig), fill=TRUE)
           cat("\tmean =", signif(x$mean, sigdig), fill=TRUE)
         },
         {
           print(x$summary)
         })

  return(invisible(NULL))
}

print.im <- function(x, ...) {
  cat(paste0(x$type, "-valued"), "pixel image", fill=TRUE)
  if(x$type == "factor") {
    cat("factor levels:\n")
    print(levels(x))
  }
  unitinfo <- summary(unitname(x))
  di <- x$dim
  cat(di[1], "x", di[2], "pixel array (ny, nx)", fill=TRUE)
  cat("enclosing rectangle:",
      prange(signif(x$xrange, 5)),
      "x",
      prange(signif(x$yrange, 5)),
      unitinfo$plural,
      unitinfo$explain,
      fill=TRUE)
  return(invisible(NULL))
}
