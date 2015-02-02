#
#    summary.im.R
#
#    summary() method for class "im"
#
#    $Revision: 1.19 $   $Date: 2015/02/01 03:59:04 $
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
  splat(paste0(x$type, "-valued"), "pixel image")
  unitinfo <- summary(x$units)
  pluralunits <- unitinfo$plural
  sigdig <- getOption('digits')
  di <- x$dim
  win <- x$window
  splat(di[1], "x", di[2], "pixel array (ny, nx)")
  splat("enclosing rectangle:",
        prange(signif(x$window$xrange, sigdig)),
        "x",
        prange(signif(x$window$yrange, sigdig)),
        unitinfo$plural,
        unitinfo$explain)
  splat("dimensions of each pixel:",
        signif(x$xstep, 3), "x", signif(x$ystep, sigdig),
        pluralunits)
  if(!is.null(explain <- unitinfo$explain))
    splat(explain)
  if(x$fullgrid) {
    splat("Image is defined on the full rectangular grid")
    whatpart <- "Frame"
  } else {
    splat("Image is defined on a subset of the rectangular grid")
    whatpart <- "Subset"
  }
  splat(whatpart, "area =", win$area, "square", pluralunits)
  if(x$fullgrid) splat("Pixel values") else
                 splat("Pixel values (inside window):")
  switch(x$type,
         integer=,
         real={
           splat("\trange =", prange(signif(x$range, sigdig)))
           splat("\tintegral =", signif(x$integral, sigdig))
           splat("\tmean =", signif(x$mean, sigdig))
         },
         factor={
           print(x$table)
         },
         complex={
           splat("\trange: Real",
                 prange(signif(x$Re$range, sigdig)),
                 "Imaginary",
                 prange(signif(x$Im$range, sigdig)))
           splat("\tintegral =", signif(x$integral, sigdig))
           splat("\tmean =", signif(x$mean, sigdig))
         },
         {
           print(x$summary)
         })

  return(invisible(NULL))
}

print.im <- function(x, ...) {
  splat(paste0(x$type, "-valued"), "pixel image")
  if(x$type == "factor") {
    splat("factor levels:")
    print(levels(x))
  }
  sigdig <- min(5, getOption('digits'))
  unitinfo <- summary(unitname(x))
  di <- x$dim
  splat(di[1], "x", di[2], "pixel array (ny, nx)")
  splat("enclosing rectangle:",
        prange(signif(zapsmall(x$xrange), sigdig)),
        "x",
        prange(signif(zapsmall(x$yrange), sigdig)),
        unitinfo$plural,
        unitinfo$explain)
  return(invisible(NULL))
}
