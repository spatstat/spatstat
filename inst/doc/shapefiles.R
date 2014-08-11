### R code from vignette source 'shapefiles.Rnw'

###################################################
### code chunk number 1: shapefiles.Rnw:7-8
###################################################
options(SweaveHooks=list(fig=function() par(mar=c(1,1,1,1))))


###################################################
### code chunk number 2: shapefiles.Rnw:23-29
###################################################
library(spatstat)
options(useFancyQuotes=FALSE)
sdate <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Date")
sversion <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Version")


###################################################
### code chunk number 3: shapefiles.Rnw:104-105 (eval = FALSE)
###################################################
## library(maptools)


###################################################
### code chunk number 4: shapefiles.Rnw:109-110 (eval = FALSE)
###################################################
## x <- readShapeSpatial("mydata.shp")


###################################################
### code chunk number 5: shapefiles.Rnw:115-116 (eval = FALSE)
###################################################
## class(x)


###################################################
### code chunk number 6: shapefiles.Rnw:133-137
###################################################
baltim <- columbus <- fylk <- list()
class(baltim) <- "SpatialPointsDataFrame"
class(columbus) <- "SpatialPolygonsDataFrame"
class(fylk) <- "SpatialLinesDataFrame"


###################################################
### code chunk number 7: shapefiles.Rnw:139-143 (eval = FALSE)
###################################################
## setwd(system.file("shapes", package="maptools"))
## baltim   <- readShapeSpatial("baltim.shp")
## columbus <- readShapeSpatial("columbus.shp")
## fylk     <- readShapeSpatial("fylk-val.shp")


###################################################
### code chunk number 8: shapefiles.Rnw:145-148
###################################################
class(baltim)
class(columbus)
class(fylk)


###################################################
### code chunk number 9: shapefiles.Rnw:170-171 (eval = FALSE)
###################################################
## X <- X[W]


###################################################
### code chunk number 10: shapefiles.Rnw:188-189 (eval = FALSE)
###################################################
## y <- as(x, "ppp")


###################################################
### code chunk number 11: shapefiles.Rnw:200-202 (eval = FALSE)
###################################################
## balt <- as(baltim, "ppp")
## bdata <- slot(baltim, "data")


###################################################
### code chunk number 12: shapefiles.Rnw:250-251 (eval = FALSE)
###################################################
## out <- lapply(x@lines, function(z) { lapply(z@Lines, as.psp) })


###################################################
### code chunk number 13: shapefiles.Rnw:260-261 (eval = FALSE)
###################################################
## curvegroup <- lapply(out, function(z) { do.call("superimposePSP", z)})


###################################################
### code chunk number 14: shapefiles.Rnw:301-305 (eval = FALSE)
###################################################
## out <- lapply(x@lines, function(z) { lapply(z@Lines, as.psp) })
## dat <- x@data
## for(i in seq(nrow(dat))) 
##   out[[i]] <- lapply(out[[i]], "marks<-", value=dat[i, , drop=FALSE])


###################################################
### code chunk number 15: shapefiles.Rnw:326-328
###################################################
getOption("SweaveHooks")[["fig"]]()
data(chorley)
plot(as.owin(chorley), lwd=3, main="polygon")


###################################################
### code chunk number 16: shapefiles.Rnw:341-343
###################################################
getOption("SweaveHooks")[["fig"]]()
data(demopat)
plot(as.owin(demopat), col="blue", main="polygonal region")


###################################################
### code chunk number 17: shapefiles.Rnw:379-382 (eval = FALSE)
###################################################
## regions <- slot(x, "polygons")
## regions <- lapply(regions, function(x) { SpatialPolygons(list(x)) })
## windows <- lapply(regions, as.owin)


###################################################
### code chunk number 18: shapefiles.Rnw:387-388 (eval = FALSE)
###################################################
## te <- tess(tiles=windows)


###################################################
### code chunk number 19: shapefiles.Rnw:424-425 (eval = FALSE)
###################################################
## y <- as(x, "SpatialPolygons")


###################################################
### code chunk number 20: shapefiles.Rnw:435-439 (eval = FALSE)
###################################################
## cp      <- as(columbus, "SpatialPolygons")
## cregions <- slot(cp, "polygons")
## cregions <- lapply(cregions, function(x) { SpatialPolygons(list(x)) })
## cwindows <- lapply(cregions, as.owin)


###################################################
### code chunk number 21: shapefiles.Rnw:449-451 (eval = FALSE)
###################################################
## ch <- hyperframe(window=cwindows)
## ch <- cbind.hyperframe(ch, columbus@data)


###################################################
### code chunk number 22: shapefiles.Rnw:471-473 (eval = FALSE)
###################################################
##   y <- as(x, "im")
##   ylist <- lapply(slot(x, "data"), function(z, y) { y[,] <- z; y }, y=y)


