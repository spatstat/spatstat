#
#    summary.im.R
#
#    summary() method for class "im"
#
#    $Revision: 1.15 $   $Date: 2011/05/18 09:15:30 $
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
  if(fak <- !is.null(lev))
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
  cat(paste(x$type, "-valued pixel image\n", sep=""))
  unitinfo <- summary(x$units)
  pluralunits <- unitinfo$plural
  di <- x$dim
  win <- x$window
  cat(paste(di[1], "x", di[2], "pixel array (ny, nx)\n"))
  cat("enclosing rectangle: ")
  cat(paste("[",
            paste(win$xrange, collapse=", "),
            "] x [",
            paste(win$yrange, collapse=", "),
            "] ",
            pluralunits, "\n", sep=""))
  cat(paste("dimensions of each pixel:",
            signif(x$xstep, 3), "x", signif(x$ystep, 3),
            pluralunits, "\n"))
  if(!is.null(explain <- unitinfo$explain))
    cat(paste(explain, "\n"))
  if(x$fullgrid) {
    cat("Image is defined on the full rectangular grid\n")
    whatpart <- "Frame"
  } else {
    cat("Image is defined on a subset of the rectangular grid\n")
    whatpart <- "Subset"
  }
  cat(paste(whatpart, "area = ", win$area, "square", pluralunits, "\n"))
  cat(paste("Pixel values ",
            if(x$fullgrid) "" else "(inside window)",
            ":\n", sep=""))
  switch(x$type,
         integer=,
         real={
           cat(paste(
                     "\trange = [",
                     paste(x$range, collapse=","),
                     "]\n",
                     "\tintegral = ",
                     x$integral,
                     "\n",
                     "\tmean = ",
                     x$mean,
                     "\n",
                     sep=""))
         },
         factor={
           print(x$table)
         },
         complex={
           cat(paste(
                     "\trange: Real [",
                     paste(x$Re$range, collapse=","),
                     "], Imaginary [",
                     paste(x$Im$range, collapse=","),
                     "]\n",
                     "\tintegral = ",
                     x$integral,
                     "\n",
                     "\tmean = ",
                     x$mean,
                     "\n",
                     sep=""))
         },
         {
           print(x$summary)
         })

  return(invisible(NULL))
}



print.im <- function(x, ...) {
  cat(paste(x$type, "-valued pixel image\n", sep=""))
  if(x$type == "factor") {
    cat("factor levels:\n")
    print(levels(x))
  }
  unitinfo <- summary(unitname(x))
  di <- x$dim
  cat(paste(di[1], "x", di[2], "pixel array (ny, nx)\n"))
  cat("enclosing rectangle: ")
  cat(paste("[",
            paste(signif(x$xrange, 5), collapse=", "),
            "] x [",
            paste(signif(x$yrange, 5), collapse=", "),
            "] ",
            unitinfo$plural,
            " ", unitinfo$explain,
            "\n", sep=""))
  return(invisible(NULL))
}
