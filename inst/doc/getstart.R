### R code from vignette source 'getstart.Rnw'

###################################################
### code chunk number 1: getstart.Rnw:5-6
###################################################
options(SweaveHooks=list(fig=function() par(mar=c(1,1,1,1))))


###################################################
### code chunk number 2: getstart.Rnw:23-30
###################################################
library(spatstat)
spatstat.options(image.colfun=function(n) { grey(seq(0,1,length=n)) })
sdate <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Date")
sversion <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Version")
options(useFancyQuotes=FALSE)


###################################################
### code chunk number 3: getstart.Rnw:54-56
###################################################
getOption("SweaveHooks")[["fig"]]()
data(redwood)
plot(redwood, pch=16, main="")


###################################################
### code chunk number 4: getstart.Rnw:77-79
###################################################
getOption("SweaveHooks")[["fig"]]()
data(longleaf)
plot(longleaf, main="")


###################################################
### code chunk number 5: getstart.Rnw:136-139
###################################################
data(finpines)
mypattern <- unmark(finpines)
mydata <- round(as.data.frame(finpines), 2)


###################################################
### code chunk number 6: getstart.Rnw:154-155 (eval = FALSE)
###################################################
## mydata <- read.csv("myfile.csv")


###################################################
### code chunk number 7: getstart.Rnw:165-166
###################################################
head(mydata)


###################################################
### code chunk number 8: getstart.Rnw:181-182 (eval = FALSE)
###################################################
##   mypattern <- ppp(mydata[,3], mydata[,7], c(100,200), c(10,90))


###################################################
### code chunk number 9: getstart.Rnw:185-186 (eval = FALSE)
###################################################
## ppp(x.coordinates, y.coordinates, x.range, y.range)


###################################################
### code chunk number 10: getstart.Rnw:195-196
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(mypattern)


###################################################
### code chunk number 11: getstart.Rnw:203-204 (eval = FALSE)
###################################################
## summary(mypattern)


###################################################
### code chunk number 12: getstart.Rnw:208-209
###################################################
options(SweaveHooks=list(fig=function() par(mar=rep(4,4)+0.1)))


###################################################
### code chunk number 13: getstart.Rnw:211-212
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(Kest(mypattern))


###################################################
### code chunk number 14: getstart.Rnw:218-219 (eval = FALSE)
###################################################
## plot(envelope(mypattern,Kest))


###################################################
### code chunk number 15: getstart.Rnw:221-222
###################################################
env <- envelope(mypattern,Kest, nsim=39)


###################################################
### code chunk number 16: getstart.Rnw:224-225
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(env, main="envelope(mypattern, Kest)")


###################################################
### code chunk number 17: getstart.Rnw:227-228
###################################################
options(SweaveHooks=list(fig=function() par(mar=c(1,1,1,1))))


###################################################
### code chunk number 18: getstart.Rnw:234-235
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(density(mypattern))


###################################################
### code chunk number 19: getstart.Rnw:245-246 (eval = FALSE)
###################################################
## marks(mypattern) <- mydata[, c(5,9)]


###################################################
### code chunk number 20: getstart.Rnw:248-249
###################################################
mypattern <-finpines


###################################################
### code chunk number 21: getstart.Rnw:252-253 (eval = FALSE)
###################################################
## plot(Smooth(mypattern))


###################################################
### code chunk number 22: getstart.Rnw:256-257
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(Smooth(mypattern, sigma=1.2), main="Smooth(mypattern)")


