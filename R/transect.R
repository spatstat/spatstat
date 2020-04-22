#
#  transect.R
#
# Line transects of pixel images
#
#  $Revision: 1.7 $  $Date: 2020/04/22 05:18:40 $
#

transect.im <- local({

  specify.location <- function(loc, rect) {
    lname <- short.deparse(substitute(loc))
    if(is.numeric(loc) && length(loc) == 2)
      return(list(x=loc[1], y=loc[2]))
    if(is.list(loc))
      return(xy.coords(loc))
    if(!(is.character(loc) && length(loc) == 1))
      stop(paste("Unrecognised format for", sQuote(lname)), call.=FALSE)
    xr <- rect$xrange
    yr <- rect$yrange
    switch(loc,
           bottomleft  = list(x=xr[1],    y=yr[1]),
           bottom      = list(x=mean(xr), y=yr[1]),
           bottomright = list(x=xr[2],    y=yr[1]),
           right       = list(x=xr[2],    y=mean(yr)),
           topright    = list(x=xr[2],    y=yr[2]),
           top         = list(x=mean(xr), y=yr[2]),
           topleft     = list(x=xr[1],    y=yr[2]),
           left        = list(x=xr[1],    y=mean(yr)),
           centre=,
           center      = list(x=mean(xr), y=mean(yr)),
           stop(paste("Unrecognised location",
                      sQuote(lname), "=", dQuote(loc)),
                call.=FALSE)
           )
  }

  transect.im <- 
    function(X, ..., from="bottomleft", to="topright",
             nsample=512, click=FALSE, add=FALSE) {
      Xname <- short.deparse(substitute(X))
      Xname <- sensiblevarname(Xname, "X")
      stopifnot(is.im(X))
      check.1.integer(nsample)
      # determine transect position
      if(click) {
        # interactive
        if(!add) plot(X)
        from <- spatstatLocator(1)
        points(from)
        to <- spatstatLocator(1)
        points(to)
        segments(from$x, from$y, to$x, to$y)
      } else {
        # data defining a line segment
        R <- as.rectangle(X)
        from <- specify.location(from, R)
        to   <- specify.location(to,   R)
      }
      # create sample points along transect
      if(identical(from,to))
        stop(paste(sQuote("from"), "and", sQuote("to"),
                   "must be distinct points"), call.=FALSE)
      u <- seq(0,1,length.out=nsample)
      x <- from$x + u * (to$x - from$x)
      y <- from$y + u * (to$y - from$y)
      leng <- sqrt( (to$x - from$x)^2 +  (to$y - from$y)^2)
      t <- u * leng
      # look up pixel values (may be NA)
      v <- X[list(x=x, y=y), drop=FALSE]
      # package into fv object
      df <- data.frame(t=t, v=v)
      colnames(df)[2] <- Xname
      fv(df, argu = "t",
         ylab = substitute(Xname(t), list(Xname=as.name(Xname))),
         valu=Xname,
         labl = c("t", "%s(t)"),
         desc = c("distance along transect",
           "pixel value of %s"),
         unitname = unitname(X), fname = Xname)
    }

  transect.im
})
