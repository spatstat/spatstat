### R code from vignette source 'updates.Rnw'

###################################################
### code chunk number 1: updates.Rnw:20-24
###################################################
library(spatstat)
sversion <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Version")
options(useFancyQuotes=FALSE)


###################################################
### code chunk number 2: updates.Rnw:43-49
###################################################
z <- read.table("packagesizes.txt", header=TRUE)
z$date <- as.Date(z$date)
Plot <- function(fmla, ..., dat=z) {
  yvals <- eval(as.expression(fmla[[2]]), envir=dat)
  plot(fmla, ..., data=dat, type="l", xlab="", lwd=2, ylim=c(0, max(yvals)))
}


###################################################
### code chunk number 3: updates.Rnw:51-52
###################################################
options(SweaveHooks=list(fig=function() par(mar=0.2+c(2,4,2,0))))


###################################################
### code chunk number 4: updates.Rnw:58-63
###################################################
getOption("SweaveHooks")[["fig"]]()
Plot((Rlines + srclines)/1000 ~ date, ylab="Lines of code (x 1000)", 
     main="Spatstat growth")
lines(srclines/1000 ~ date, data=z)
text(as.Date("2013-01-01"), 9.5, "C code")
text(as.Date("2013-01-01"), 50, "R code")


###################################################
### code chunk number 5: updates.Rnw:2616-2620
###################################################
nbugs <- nrow(news(grepl("^BUG", Category), 
                   package="spatstat"))
nbugssince <- nrow(news(Version > "1.21-2" & grepl("^BUG", Category), 
                   package="spatstat"))


###################################################
### code chunk number 6: updates.Rnw:2627-2628 (eval = FALSE)
###################################################
## news(grepl("^BUG", Category), package="spatstat")


