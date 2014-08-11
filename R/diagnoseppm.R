#
#                            diagnoseppm.R
#
# Makes diagnostic plots based on residuals or energy weights
#
# $Revision: 1.32 $ $Date: 2011/02/15 07:41:16 $
#

diagnose.ppm.engine <- function(object, ..., type="eem", typename, opt,
                              sigma=NULL,
                              rbord = reach(object), compute.sd=TRUE,
                              compute.cts=TRUE, rv=NULL, oldstyle=FALSE)
{
  if(is.marked.ppm(object))
    stop("Sorry, this is not yet implemented for marked models")

  # quadrature points
  Q <- quad.ppm(object)
  U <- union.quad(Q)
  Qweights <- w.quad(Q)
  
  # -------------- Calculate residuals/weights -------------------

  # Discretised residuals

  if(type == "eem") {
    residval <- if(!is.null(rv)) rv else eem(object, check=FALSE)
    residval <- as.numeric(residval)
    X <- data.ppm(object)
    Y <- X %mark% residval
  } else {
    if(!is.null(rv) && !inherits(rv, "msr"))
      stop("rv should be a measure (object of class msr)")
    residobj <-
      if(!is.null(rv)) rv else residuals.ppm(object, type=type, check=FALSE)
    residval <- with(residobj, "increment")
    if(ncol(as.matrix(residval)) > 1)
      stop("Not implemented for vector-valued residuals; use [.msr to split into separate components")
    Y <- U %mark% residval
  }

  # Atoms and density of measure

  Ymass <- NULL
  Ycts  <- NULL
  Ydens <- NULL

  if(compute.cts) {
    if(type == "eem") {
      Ymass <- Y
      Ycts  <- U %mark% (-1)
      Ydens <- as.im(-1, Y$window)
    } else {
      atoms <- with(residobj, "is.atom")
      masses <- with(residobj, "discrete")
      cts    <- with(residobj, "density")
      if(!is.null(atoms) && !is.null(masses) && !is.null(cts)) {
        Ymass <- (U %mark% masses)[atoms]
        Ycts    <- U %mark% cts
        # remove NAs (as opposed to zero cif points)
        if(!all(ok <- is.finite(cts))) {
          U <- U[ok]
          Ycts <- Ycts[ok]
          cts  <- cts[ok]
          Qweights <- Qweights[ok]
        }
        # interpolate continuous part to yield an image for plotting
        if(type == "inverse" && all(cts > 0)) {
          Ydens <- as.im(-1, Y$window)
        } else if(is.stationary.ppm(object) && is.poisson.ppm(object)) {
          # all values of `cts' will be equal
          Ydens <- as.im(cts[1], Y$window)
        } else {
          smallsigma <- max(nndist(Ycts))
          Ujitter <- U
          Ujitter$x <- U$x + runif(U$n, -smallsigma, smallsigma)
          Ujitter$y <- U$y + runif(U$n, -smallsigma, smallsigma)
          Ydens <- smooth.ppp(Ujitter %mark% marks(Ycts),
                              sigma=smallsigma,
                              weights=Qweights,
                              edge=TRUE, ...)
        }
      }
    }
  }
    

  #----------------  Erode window ---------------------------------
  #
  ## Compute windows 
  W <- Y$window

  # Erode window if required
  clip <- (rbord > 0)
  if(clip) {
    Wclip <- erosion.owin(W, rbord)
    Yclip <- Y[Wclip]
    Qweightsclip <- Qweights[inside.owin(U, , Wclip)]
    if(!is.null(Ycts))
      Ycts <- Ycts[Wclip]
    if(!is.null(Ydens))
      Ydens <- Ydens[Wclip, drop=FALSE]
  } else {
    Wclip <- W
    Yclip <- Y
  }
  
  # ------------ start collecting results -------------------------
  
  result <- list(type=type,
                 clip=clip,
                 Y=Y,
                 W=W,
                 Yclip=Yclip,
                 Ymass=Ymass,
                 Ycts=Ycts,
                 Ydens=Ydens)

  # ------------- smoothed field ------------------------------

  Z <- NULL
  if(opt$smooth | opt$xcumul | opt$ycumul | opt$xmargin | opt$ymargin) {
    if(is.null(sigma))
      sigma <- 0.1 * diameter(Wclip)  
    Z <- density.ppp(Yclip, sigma, weights=Yclip$marks, edge=TRUE, ...)
  }
  if(opt$smooth)
    result$smooth <- list(Z = Z, sigma=sigma)

  # -------------- marginals of smoothed field ------------------------
  
  if(opt$xmargin) {
    xZ <- apply(Z$v, 2, sum, na.rm=TRUE) * Z$xstep
    if(type == "eem") 
      ExZ <- apply(Z$v, 2, function(column) { sum(!is.na(column)) }) * Z$xstep
    else 
      ExZ <- rep(0, length(xZ))
    result$xmargin <- list(x=Z$xcol, xZ=xZ, ExZ=ExZ)
  }
  
  if(opt$ymargin) {
    yZ <- apply(Z$v, 1, sum, na.rm=TRUE) * Z$ystep
    if(type == "eem")
      EyZ <- apply(Z$v, 1, function(roww) { sum(!is.na(roww)) }) * Z$ystep
    else
      EyZ <- rep(0, length(yZ))
    result$ymargin <- list(y=Z$yrow, yZ=yZ, EyZ=EyZ)
  }
  
  # -------------- cumulative (lurking variable) plots --------------

  if(opt$xcumul)
    result$xcumul <- 
    lurking(object, covariate=x.quad(Q),
            type=type,
            clipwindow= if(clip) Wclip else NULL,
            rv=residval,
            plot.sd=compute.sd,
            plot.it=FALSE,
            typename=typename,
            covname="x coordinate",
            oldstyle=oldstyle,
            check=FALSE, ...)

  if(opt$ycumul)
    result$ycumul <- 
    lurking(object, covariate=y.quad(Q),
            type=type,
            clipwindow= if(clip) Wclip else NULL,
            rv=residval,
            plot.sd=compute.sd,
            plot.it=FALSE,
            typename=typename,
            covname="y coordinate",
            oldstyle=oldstyle,
            check=FALSE, ...)

  # -------------- summary numbers --------------
  
  if(opt$sum) 
    result$sum <- list(marksum=sum(Yclip$marks, na.rm=TRUE),
                       area.Wclip=area.owin(Wclip),
                       area.quad=if(clip) sum(Qweightsclip) else sum(Qweights),
                       range=if(!is.null(Z)) range(Z) else NULL)

  return(invisible(result))
}


########################################################################


diagnose.ppm <- function(object, ..., type="raw", which="all", 
                         sigma=NULL, 
                         rbord = reach(object), cumulative=TRUE,
                         plot.it = TRUE, rv = NULL, 
                         compute.sd=TRUE, compute.cts=TRUE,
                         typename, check=TRUE, repair=TRUE, oldstyle=FALSE)
{
  if(is.marked.ppm(object))
    stop("Sorry, this is not yet implemented for marked models")

  if(check && damaged.ppm(object)) {
    if(!repair)
      stop("object format corrupted; try update(object, use.internal=TRUE)")
    message("object format corrupted; repairing it.")
    object <- update(object, use.internal=TRUE)
  }
    
  # -------------  Interpret arguments --------------------------

  # edge effect avoidance
  if(!is.finite(rbord)) {
    if(missing(rbord))
      stop(paste(sQuote("rbord"),
                 "must be specified; the model has infinite range"))
    else
      stop(paste(sQuote("rbord"), "is infinite"))
  }
  
  # whether window should be clipped
  clip <- (rbord > 0)

  # match type argument
  type <- pickoption("type", type,
                     c(eem="eem",
                       raw="raw",
                       inverse="inverse",
                       pearson="pearson",
                       Pearson="pearson"))
  if(missing(typename))
    typename <- switch(type,
                       eem="exponential energy weights",
                       raw="raw residuals",
                       inverse="inverse-lambda residuals",
                       pearson="Pearson residuals")

  # 'which' is multiple choice with exact matching 
  optionlist <- c("all", "marks", "smooth", "x", "y", "sum")

  if(!all(m <- which %in% optionlist))
    stop(paste("Unrecognised choice(s) of",
               paste(sQuote("which"), ":", sep=""),
               paste(which[!m], collapse=", ")))

  opt <- list()
  opt$all <- "all" %in% which
  opt$marks <-  ("marks" %in% which)   | opt$all
  opt$smooth <- ("smooth" %in% which)  | opt$all
  opt$xmargin <- (("x" %in% which)       | opt$all) && !cumulative
  opt$ymargin <- (("y" %in% which)       | opt$all) && !cumulative
  opt$xcumul <-  (("x" %in% which)       | opt$all) && cumulative
  opt$ycumul <-  (("y" %in% which)       | opt$all) && cumulative
  opt$sum <-     ("sum" %in% which)      | opt$all

  # compute and plot estimated standard deviations?
  # yes for Poisson, no for other models, unless overridden
  if(!missing(compute.sd))
    plot.sd <- compute.sd
  else
    plot.sd <- list(...)$plot.sd
  if(is.null(plot.sd))
    plot.sd <- is.poisson.ppm(object)
  if(missing(compute.sd))
    compute.sd <- plot.sd

  # interpolate the density of the residual measure?
  if(missing(compute.cts)) {
    plot.neg <- resolve.defaults(list(...),
                                 formals(plot.diagppm)["plot.neg"])$plot.neg
    # only if it is needed for the mark plot
    compute.cts <- opt$marks && (plot.neg != "discrete")
  }

  # -------  DO THE CALCULATIONS -----------------------------------
  RES <-  diagnose.ppm.engine(object, type=type, typename=typename,
                              opt=opt, sigma=sigma, rbord=rbord,
                              compute.sd=compute.sd,
                              compute.cts=compute.cts,
                              rv=rv, oldstyle=oldstyle, ...)

  RES$typename <- typename
  RES$opt <- opt
  RES$compute.sd <- compute.sd
  RES$compute.cts <- compute.cts
  
  class(RES) <- "diagppm"

  # -------  PLOT --------------------------------------------------
  if(plot.it)
    plot(RES, ...)

  return(RES)
}

plot.diagppm <- function(x, ..., which, plot.neg="image",
                         plot.smooth="imagecontour",
                         plot.sd=TRUE, spacing=0.1,
                         srange=NULL, monochrome=FALSE, main=NULL)
{
  opt <- x$opt

  if(!missing(which)) {
    oldopt <- opt
    newopt <- list()
    newopt$all <- "all" %in% which
    newopt$marks <-  ("marks" %in% which)   | newopt$all
    newopt$smooth <- ("smooth" %in% which)  | newopt$all
    newopt$xmargin <- (("x" %in% which)       | newopt$all) && oldopt$xmargin
    newopt$ymargin <- (("y" %in% which)       | newopt$all) && oldopt$ymargin
    newopt$xcumul <-  (("x" %in% which)       | newopt$all) && oldopt$xcumul
    newopt$ycumul <-  (("y" %in% which)       | newopt$all)  && oldopt$ycumul
    newopt$sum <-     ("sum" %in% which)      | newopt$all

    illegal <- (unlist(newopt) > unlist(oldopt))
    if(any(illegal)) {
      offending <- paste(names(newopt)[illegal], collapse=", ")
      whinge <- paste("cannot display the following components;\n",
                      "they were not computed: - \n", offending, "\n")
      stop(whinge)
    }

    opt <- newopt
  }

  if(!(x$compute.sd) && plot.sd) {
    if(!missing(plot.sd))
      warning("can't plot standard deviations; they were not computed")
    plot.sd <- FALSE
  }

  if(!(x$compute.cts) && (plot.neg != "discrete") && (opt$marks || opt$all)) {
    if(!missing(plot.neg))
      warning("can't plot continuous component of residuals; it was not computed")
    plot.neg <- "discrete"
  }
  
  if(opt$all) 
    resid4plot(x, plot.neg, plot.smooth, spacing, srange,monochrome, main, ...)
  else
    resid1plot(x, opt, plot.neg, plot.smooth, srange, monochrome, main, ...)
}


print.diagppm <- function(x, ...) {
  
  opt <- x$opt
  typename <- x$typename
  
  cat(paste("Model diagnostics (", typename, ")\n", sep=""))

  cat("Diagnostics available:\n")
  optkey <- list(all="four-panel plot",
                 marks=paste("mark plot", if(!x$compute.cts)
                   "(discrete representation only)" else NULL),
                 smooth="smoothed residual field",
                 xmargin="x marginal density",
                 ymargin="y marginal density",
                 xcumul="x cumulative residuals",
                 ycumul="y cumulative residuals",
                 sum="sum of all residuals")
  avail <- unlist(optkey[names(opt)[unlist(opt)]])
  names(avail) <- NULL
  cat(paste("\t", paste(avail, collapse="\n\t"), "\n", sep=""))
  
  if(opt$sum) {
    xs <- x$sum
    windowname <- if(x$clip) "clipped window" else "entire window"
    cat(paste("sum of", typename, "in", windowname, "=",
              signif(sum(xs$marksum),4), "\n"))
    cat(paste("area of", windowname, "=",
              signif(xs$area.Wclip, 4), "\n"))
    cat(paste("quadrature area =",
              signif(xs$area.quad, 4), "\n"))
  }
  if(opt$smooth) 
    cat(paste("range of smoothed field = [",
              paste(signif(range(x$smooth$Z$v, na.rm=TRUE),4), collapse=","),
              "]\n"))

  return(invisible(NULL))
}
