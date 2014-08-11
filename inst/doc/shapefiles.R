### R code from vignette source 'shapefiles.Rnw'

###################################################
### code chunk number 1: shapefiles.Rnw:7-8
###################################################
options(SweaveHooks=list(fig=function() par(mar=c(1,1,1,1))))


###################################################
### code chunk number 2: shapefiles.Rnw:23-30
###################################################
library(spatstat)
spatstat.options(gpclib=TRUE)
options(useFancyQuotes=FALSE)
sdate <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Date")
sversion <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Version")


###################################################
### code chunk number 3: shapefiles.Rnw:113-114 (eval = FALSE)
###################################################
## library(maptools)


###################################################
### code chunk number 4: shapefiles.Rnw:118-119 (eval = FALSE)
###################################################
## x <- readShapeSpatial("mydata.shp")


###################################################
### code chunk number 5: shapefiles.Rnw:124-125 (eval = FALSE)
###################################################
## class(x)


###################################################
### code chunk number 6: shapefiles.Rnw:140-144
###################################################
baltim <- columbus <- fylk <- list()
class(baltim) <- "SpatialPointsDataFrame"
class(columbus) <- "SpatialPolygonsDataFrame"
class(fylk) <- "SpatialLinesDataFrame"


###################################################
### code chunk number 7: shapefiles.Rnw:146-150 (eval = FALSE)
###################################################
## setwd(system.file("shapes", package="maptools"))
## baltim   <- readShapeSpatial("baltim.shp")
## columbus <- readShapeSpatial("columbus.shp")
## fylk     <- readShapeSpatial("fylk-val.shp")


###################################################
### code chunk number 8: shapefiles.Rnw:152-155
###################################################
class(baltim)
class(columbus)
class(fylk)


###################################################
### code chunk number 9: shapefiles.Rnw:177-178 (eval = FALSE)
###################################################
## X <- X[W]


###################################################
### code chunk number 10: shapefiles.Rnw:195-196 (eval = FALSE)
###################################################
## y <- as(x, "ppp")


###################################################
### code chunk number 11: shapefiles.Rnw:207-209 (eval = FALSE)
###################################################
## balt <- as(baltim, "ppp")
## bdata <- slot(baltim, "data")


###################################################
### code chunk number 12: shapefiles.Rnw:257-258 (eval = FALSE)
###################################################
## out <- lapply(x@lines, function(z) { lapply(z@Lines, as.psp) })


###################################################
### code chunk number 13: shapefiles.Rnw:267-268 (eval = FALSE)
###################################################
## curvegroup <- lapply(out, function(z) { do.call("superimposePSP", z)})


###################################################
### code chunk number 14: shapefiles.Rnw:308-312 (eval = FALSE)
###################################################
## out <- lapply(x@lines, function(z) { lapply(z@Lines, as.psp) })
## dat <- x@data
## for(i in seq(nrow(dat))) 
##   out[[i]] <- lapply(out[[i]], "marks<-", value=dat[i, , drop=FALSE])


###################################################
### code chunk number 15: shapefiles.Rnw:333-335
###################################################
getOption("SweaveHooks")[["fig"]]()
data(chorley)
plot(as.owin(chorley), lwd=3, main="polygon")


###################################################
### code chunk number 16: shapefiles.Rnw:348-350
###################################################
getOption("SweaveHooks")[["fig"]]()
data(demopat)
plot(as.owin(demopat), col="blue", main="polygonal region")


###################################################
### code chunk number 17: shapefiles.Rnw:386-389 (eval = FALSE)
###################################################
## regions <- slot(x, "polygons")
## regions <- lapply(regions, function(x) { SpatialPolygons(list(x)) })
## windows <- lapply(regions, as.owin)


###################################################
### code chunk number 18: shapefiles.Rnw:394-395 (eval = FALSE)
###################################################
## te <- tess(tiles=windows)


###################################################
### code chunk number 19: shapefiles.Rnw:427-430 (eval = FALSE)
###################################################
## spatstat.options(checkpolygons=FALSE)
## y <- as(x, "owin")
## spatstat.options(checkpolygons=TRUE)


###################################################
### code chunk number 20: shapefiles.Rnw:447-448 (eval = FALSE)
###################################################
## y <- as(x, "SpatialPolygons")


###################################################
### code chunk number 21: shapefiles.Rnw:458-462 (eval = FALSE)
###################################################
## cp      <- as(columbus, "SpatialPolygons")
## cregions <- slot(cp, "polygons")
## cregions <- lapply(cregions, function(x) { SpatialPolygons(list(x)) })
## cwindows <- lapply(cregions, as.owin)


###################################################
### code chunk number 22: shapefiles.Rnw:472-474 (eval = FALSE)
###################################################
## ch <- hyperframe(window=cwindows)
## ch <- cbind.hyperframe(ch, columbus@data)


