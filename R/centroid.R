#
#	centroid.S	Centroid of a window
#			and related operations
#
#	$Revision: 1.2 $	$Date: 2013/05/01 05:42:54 $
#
# Function names (followed by "xypolygon" or "owin")
#	
#	intX            integral of x dx dy
#	intY            integral of y dx dy
#	meanX           mean of x dx dy
#	meanY           mean of y dx dy
#       centroid        (meanX, meanY)
#		
#-------------------------------------

intX.xypolygon <- function(polly) {
  #
  # polly: list(x,y) vertices of a single polygon (n joins to 1)
  #
  verify.xypolygon(polly)
  
  x <- polly$x
  y <- polly$y
  
  nedges <- length(x)   # sic
  
  # place x axis below polygon
  y <- y - min(y) 

  # join vertex n to vertex 1
  xr <- c(x, x[1])
  yr <- c(y, y[1])

  # slope
  dx <- diff(xr)
  dy <- diff(yr)
  slope <- ifelseAX(dx == 0, 0, dy/dx)

  # integrate
  integrals <- x * y * dx + (y + slope * x) * (dx^2)/2 + slope * (dx^3)/3

  -sum(integrals)
}
		
intX.owin <- function(w) {
	verifyclass(w, "owin")
        switch(w$type,
               rectangle = {
		width  <- abs(diff(w$xrange))
		height <- abs(diff(w$yrange))
		answer <- width * height * mean(w$xrange)
               },
               polygonal = {
                 answer <- sum(unlist(lapply(w$bdry, intX.xypolygon)))
               },
               mask = {
                 pixelarea <- abs(w$xstep * w$ystep)
		 npixels <- sum(w$m)
		 area <- npixels * pixelarea
		 x <- raster.x(w)[w$m]
                 answer <- area * mean(x)
               },
               stop("Unrecognised window type")
        )
        return(answer)
}

meanX.owin <- function(w) {
	verifyclass(w, "owin")
        switch(w$type,
               rectangle = {
		answer <- mean(w$xrange)
               },
               polygonal = {
	         area <- sum(unlist(lapply(w$bdry, area.xypolygon)))
                 integrated <- sum(unlist(lapply(w$bdry, intX.xypolygon)))
		 answer <- integrated/area
               },
               mask = {
		 x <- raster.x(w)[w$m]
                 answer <- mean(x)
               },
               stop("Unrecognised window type")
        )
        return(answer)
}

intY.xypolygon <- function(polly) {
  #
  # polly: list(x,y) vertices of a single polygon (n joins to 1)
  #
  verify.xypolygon(polly)
  
  x <- polly$x
  y <- polly$y
  
  nedges <- length(x)   # sic
  
  # place x axis below polygon
  yadjust <- min(y)
  y <- y - yadjust 

  # join vertex n to vertex 1
  xr <- c(x, x[1])
  yr <- c(y, y[1])

  # slope
  dx <- diff(xr)
  dy <- diff(yr)
  slope <- ifelseAX(dx == 0, 0, dy/dx)

  # integrate
  integrals <- (1/2) * (dx * y^2 + slope * y * dx^2 + slope^2 * dx^3/3)
  total <- sum(integrals) - yadjust * area.xypolygon(polly)

  # change sign to adhere to anticlockwise convention
  -total
}
		
intY.owin <- function(w) {
	verifyclass(w, "owin")
        switch(w$type,
               rectangle = {
		width  <- abs(diff(w$xrange))
		height <- abs(diff(w$yrange))
		answer <- width * height * mean(w$yrange)
               },
               polygonal = {
                 answer <- sum(unlist(lapply(w$bdry, intY.xypolygon)))
               },
               mask = {
                 pixelarea <- abs(w$xstep * w$ystep)
		 npixels <- sum(w$m)
		 area <- npixels * pixelarea
		 y <- raster.y(w)[w$m]
                 answer <- area * mean(y)
               },
               stop("Unrecognised window type")
        )
        return(answer)
}

meanY.owin <- function(w) {
	verifyclass(w, "owin")
        switch(w$type,
               rectangle = {
		answer <- mean(w$yrange)
               },
               polygonal = {
	         area <- sum(unlist(lapply(w$bdry, area.xypolygon)))
                 integrated <- sum(unlist(lapply(w$bdry, intY.xypolygon)))
		 answer <- integrated/area
               },
               mask = {
		 y <- raster.y(w)[w$m]
                 answer <- mean(y)
               },
               stop("Unrecognised window type")
        )
        return(answer)
}

centroid.owin <- function(w) {
	verifyclass(w, "owin")
	return(list(x=meanX.owin(w), y=meanY.owin(w)))
}

	
