### R code from vignette source 'replicated.Rnw'

###################################################
### code chunk number 1: replicated.Rnw:29-30
###################################################
options(SweaveHooks=list(fig=function() par(mar=c(1,1,1,1))))


###################################################
### code chunk number 2: replicated.Rnw:35-42
###################################################
library(spatstat)
spatstat.options(image.colfun=function(n) { grey(seq(0,1,length=n)) })
sdate <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Date")
sversion <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Version")
options(useFancyQuotes=FALSE)


###################################################
### code chunk number 3: replicated.Rnw:189-190
###################################################
waterstriders


###################################################
### code chunk number 4: replicated.Rnw:208-209
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(waterstriders, main="")


###################################################
### code chunk number 5: replicated.Rnw:216-217
###################################################
summary(waterstriders)


###################################################
### code chunk number 6: replicated.Rnw:225-226
###################################################
X <- listof(rpoispp(100), rpoispp(100), rpoispp(100))


###################################################
### code chunk number 7: replicated.Rnw:231-233
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(X)
X


###################################################
### code chunk number 8: replicated.Rnw:262-263 (eval = FALSE)
###################################################
## hyperframe(...)


###################################################
### code chunk number 9: replicated.Rnw:288-290
###################################################
H <- hyperframe(X=1:3, Y=list(sin,cos,tan))
H


###################################################
### code chunk number 10: replicated.Rnw:298-303
###################################################
G <- hyperframe(X=1:3, Y=letters[1:3], Z=factor(letters[1:3]),
                W=list(rpoispp(100),rpoispp(100), rpoispp(100)),
                U=42,
                V=rpoispp(100), stringsAsFactors=FALSE)
G


###################################################
### code chunk number 11: replicated.Rnw:332-333
###################################################
simba


###################################################
### code chunk number 12: replicated.Rnw:346-347
###################################################
pyramidal


###################################################
### code chunk number 13: replicated.Rnw:353-354
###################################################
ws <- hyperframe(Striders=waterstriders)


###################################################
### code chunk number 14: replicated.Rnw:361-363
###################################################
H$X
H$Y


###################################################
### code chunk number 15: replicated.Rnw:373-375
###################################################
H$U <- letters[1:3]
H


###################################################
### code chunk number 16: replicated.Rnw:380-384
###################################################
G <- hyperframe()
G$X <- waterstriders
G$Y <- 1:3
G


###################################################
### code chunk number 17: replicated.Rnw:392-396
###################################################
H[,1]
H[2,]
H[2:3, ]
H[1,1]


###################################################
### code chunk number 18: replicated.Rnw:402-405
###################################################
H[,1,drop=TRUE]
H[1,1,drop=TRUE]
H[1,2,drop=TRUE]


###################################################
### code chunk number 19: replicated.Rnw:418-419 (eval = FALSE)
###################################################
## plot.listof(x, ..., main, arrange = TRUE, nrows = NULL, ncols = NULL)


###################################################
### code chunk number 20: replicated.Rnw:434-435
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(waterstriders, pch=16, nrows=1)


###################################################
### code chunk number 21: replicated.Rnw:450-451
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(simba)


###################################################
### code chunk number 22: replicated.Rnw:463-465
###################################################
getOption("SweaveHooks")[["fig"]]()
H <- hyperframe(X=1:3, Y=list(sin,cos,tan))
plot(H$Y)


###################################################
### code chunk number 23: replicated.Rnw:477-478 (eval = FALSE)
###################################################
## plot(h, e)


###################################################
### code chunk number 24: replicated.Rnw:487-488
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(demohyper, quote({ plot(Image, main=""); plot(Points, add=TRUE) }))


###################################################
### code chunk number 25: replicated.Rnw:500-502
###################################################
getOption("SweaveHooks")[["fig"]]()
H <- hyperframe(Bugs=waterstriders)
plot(H, quote(plot(Kest(Bugs))), marsize=1)


###################################################
### code chunk number 26: replicated.Rnw:515-517
###################################################
df <- data.frame(A=1:10, B=10:1)
with(df, A-B)


###################################################
### code chunk number 27: replicated.Rnw:530-531 (eval = FALSE)
###################################################
## with(h,e)


###################################################
### code chunk number 28: replicated.Rnw:541-544
###################################################
H <- hyperframe(Bugs=waterstriders)
with(H, npoints(Bugs))
with(H, distmap(Bugs))


###################################################
### code chunk number 29: replicated.Rnw:567-568
###################################################
with(simba, npoints(Points))


###################################################
### code chunk number 30: replicated.Rnw:575-577
###################################################
H <- hyperframe(Bugs=waterstriders)
K <- with(H, Kest(Bugs))


###################################################
### code chunk number 31: replicated.Rnw:585-586
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(K)


###################################################
### code chunk number 32: replicated.Rnw:591-593
###################################################
H <- hyperframe(Bugs=waterstriders)
with(H, nndist(Bugs))


###################################################
### code chunk number 33: replicated.Rnw:599-600
###################################################
with(H, min(nndist(Bugs)))


###################################################
### code chunk number 34: replicated.Rnw:612-613
###################################################
simba$Dist <- with(simba, distmap(Points))


###################################################
### code chunk number 35: replicated.Rnw:626-630
###################################################
getOption("SweaveHooks")[["fig"]]()
lambda <- rexp(6, rate=1/50)
H <- hyperframe(lambda=lambda)
H$Points <- with(H, rpoispp(lambda))
plot(H, quote(plot(Points, main=paste("lambda=", signif(lambda, 4)))))


###################################################
### code chunk number 36: replicated.Rnw:636-637
###################################################
H$X <- with(H, rpoispp(50))


###################################################
### code chunk number 37: replicated.Rnw:666-667
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(simba, quote(plot(density(Points), main="")), nrows=2)


###################################################
### code chunk number 38: replicated.Rnw:686-688
###################################################
getOption("SweaveHooks")[["fig"]]()
rhos <- with(demohyper, rhohat(Points, Image))
plot(rhos)


###################################################
### code chunk number 39: replicated.Rnw:705-706 (eval = FALSE)
###################################################
## mppm(formula, data, interaction, ...)


###################################################
### code chunk number 40: replicated.Rnw:716-717 (eval = FALSE)
###################################################
## mppm(Points ~ group, simba, Poisson())


###################################################
### code chunk number 41: replicated.Rnw:750-751
###################################################
mppm(Points ~ 1, simba)


###################################################
### code chunk number 42: replicated.Rnw:758-759
###################################################
mppm(Points ~ group, simba)


###################################################
### code chunk number 43: replicated.Rnw:765-766
###################################################
mppm(Points ~ id, simba)


###################################################
### code chunk number 44: replicated.Rnw:776-777
###################################################
mppm(Points ~ Image, data=demohyper)


###################################################
### code chunk number 45: replicated.Rnw:795-796 (eval = FALSE)
###################################################
## mppm(Points ~ offset(log(Image)), data=demohyper)


###################################################
### code chunk number 46: replicated.Rnw:808-809 (eval = FALSE)
###################################################
## mppm(Points ~ log(Image), data=demop)


###################################################
### code chunk number 47: replicated.Rnw:826-827 (eval = FALSE)
###################################################
## mppm(formula, data, interaction, ..., iformula=NULL)


###################################################
### code chunk number 48: replicated.Rnw:877-878
###################################################
radii <- with(simba, mean(nndist(Points)))


###################################################
### code chunk number 49: replicated.Rnw:885-887
###################################################
Rad <- hyperframe(R=radii)
Str <- with(Rad, Strauss(R))


###################################################
### code chunk number 50: replicated.Rnw:892-894
###################################################
Int <- hyperframe(str=Str)
mppm(Points ~ 1, simba, interaction=Int)


###################################################
### code chunk number 51: replicated.Rnw:921-924
###################################################
h <- hyperframe(Y=waterstriders)
g <- hyperframe(po=Poisson(), str4 = Strauss(4), str7= Strauss(7))
mppm(Y ~ 1, data=h, interaction=g, iformula=~str4)


###################################################
### code chunk number 52: replicated.Rnw:935-936
###################################################
fit <- mppm(Points ~ 1, simba, Strauss(0.07), iformula = ~Interaction*group)


###################################################
### code chunk number 53: replicated.Rnw:954-955
###################################################
fit


###################################################
### code chunk number 54: replicated.Rnw:958-960
###################################################
co <- coef(fit)
si <- function(x) { signif(x, 4) }


###################################################
### code chunk number 55: replicated.Rnw:971-972
###################################################
coef(fit)


###################################################
### code chunk number 56: replicated.Rnw:1029-1030 (eval = FALSE)
###################################################
## interaction=hyperframe(po=Poisson(), str=Strauss(0.07))


###################################################
### code chunk number 57: replicated.Rnw:1035-1036 (eval = FALSE)
###################################################
## iformula=~ifelse(group=="control", po, str)


###################################################
### code chunk number 58: replicated.Rnw:1046-1047 (eval = FALSE)
###################################################
## iformula=~I((group=="control")*po) + I((group=="treatment") * str)


###################################################
### code chunk number 59: replicated.Rnw:1057-1062
###################################################
g <- hyperframe(po=Poisson(), str=Strauss(0.07))
fit2 <- mppm(Points ~ 1, simba, g, 
             iformula=~I((group=="control")*po) 
                     + I((group=="treatment") * str))
fit2


###################################################
### code chunk number 60: replicated.Rnw:1085-1088
###################################################
H <- hyperframe(W=waterstriders)
fit <- mppm(W ~ 1, H)
subfits(fit)


###################################################
### code chunk number 61: replicated.Rnw:1109-1110 (eval = FALSE)
###################################################
## subfits <- subfits.new


###################################################
### code chunk number 62: replicated.Rnw:1122-1124
###################################################
H <- hyperframe(W=waterstriders)
with(H, ppm(W))


###################################################
### code chunk number 63: replicated.Rnw:1147-1149
###################################################
fit <- mppm(P ~ x, hyperframe(P=waterstriders))
res <- residuals(fit)


###################################################
### code chunk number 64: replicated.Rnw:1159-1160
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(res)


###################################################
### code chunk number 65: replicated.Rnw:1165-1167
###################################################
getOption("SweaveHooks")[["fig"]]()
smor <- with(hyperframe(res=res), Smooth(res, sigma=4))
plot(smor)


###################################################
### code chunk number 66: replicated.Rnw:1179-1182
###################################################
fit <- mppm(P ~ x, hyperframe(P=waterstriders))
res <- residuals(fit)
totres <- sapply(res, integral.msr)


###################################################
### code chunk number 67: replicated.Rnw:1188-1195
###################################################
getOption("SweaveHooks")[["fig"]]()
fit <- mppm(Points~Image, data=demohyper)
resids <- residuals(fit, type="Pearson")
totres <- sapply(resids, integral.msr)
areas <- with(demohyper, area.owin(as.owin(Points)))
df <- as.data.frame(demohyper[, "Group"])
df$resids <- totres/areas
plot(resids~Group, df)


###################################################
### code chunk number 68: replicated.Rnw:1216-1219
###################################################
getOption("SweaveHooks")[["fig"]]()
fit <- mppm(P ~ 1, hyperframe(P=waterstriders))
sub <- hyperframe(Model=subfits(fit))
plot(sub, quote(diagnose.ppm(Model)))


###################################################
### code chunk number 69: replicated.Rnw:1232-1240
###################################################
H <- hyperframe(P = waterstriders)
fitall <- mppm(P ~ 1, H)
together <- subfits(fitall)
separate <- with(H, ppm(P))
Fits <- hyperframe(Together=together, Separate=separate)
dr <- with(Fits, unlist(coef(Separate)) - unlist(coef(Together)))
dr
exp(dr)


###################################################
### code chunk number 70: replicated.Rnw:1257-1266
###################################################
H <- hyperframe(X=waterstriders)

# Poisson with constant intensity for all patterns
fit1 <- mppm(X~1, H)
quadrat.test(fit1, nx=2)

# uniform Poisson with different intensity for each pattern
fit2 <- mppm(X ~ id, H)
quadrat.test(fit2, nx=2)


###################################################
### code chunk number 71: replicated.Rnw:1295-1296 (eval = FALSE)
###################################################
## kstest.mppm(model, covariate)


