## demonstration of all summary functions

opa <- par(mfrow=c(1,1),
           mar=c(0,0,1,0)+0.2)

## Ripley's K-function
plot(swedishpines)

par(mar=c(4,4,2,1)+0.2)
plot(Kest(swedishpines))

## Besag's transformation
plot(Lest(swedishpines))

## pair correlation function
plot(pcf(swedishpines))

par(mfrow=c(2,3),
    mar=c(0,0,1,0)+0.2)
    
## Showing the utility of the K-function
plot(cells)
plot(nztrees)
plot(redwood)
par(mar=c(4,4,2,1)+0.2)
plot(Kest(cells))
plot(Kest(nztrees))
plot(Kest(redwood))
## Showing the utility of the pair correlation function
par(mar=c(0,0,1,0)+0.2)
plot(cells)
plot(nztrees)
plot(redwood)
par(mar=c(4,4,2,1)+0.2)
plot(pcf(cells))
plot(pcf(nztrees))
plot(pcf(redwood))
##
par(mfrow=c(1,1))

## Analogues for inhomogeneous patterns
## Reweighted K-function
plot(japanesepines)
fit <- ppm(japanesepines, ~polynom(x,y,2))
plot(predict(fit))
plot(Kinhom(japanesepines, fit))
plot(pcfinhom(japanesepines, fit))
plot(Linhom(japanesepines))

## Rescaled K-function
plot(unmark(bronzefilter))
plot(Kscaled(bronzefilter))
fit <- ppm(unmark(bronzefilter), ~x)
plot(predict(fit))
plot(unmark(bronzefilter), add=TRUE)
plot(Kscaled(bronzefilter, fit))
plot(Lscaled(bronzefilter, fit))

## Local indicators of spatial association
plot(localL(swedishpines))
plot(localK(swedishpines))

## anisotropic
plot(Ksector(redwood, 0, 90))
plot(Rf <- pairorient(redwood, 0.05, 0.15))
rose(Rf, main="Rose diagram of pair orientation distribution")
plot(deriv(Rf, spar=0.6, Dperiodic=TRUE))
rose(nnorient(redwood))

##
par(mfrow=c(2,3),
    mar=rep(0.2, 4))
## Empty space function F
plot(cells)
plot(nztrees)
plot(redwood)
par(mar=c(4,4,2,1)+0.2)
plot(Fest(cells))
plot(Fest(nztrees))
plot(Fest(redwood))
## Nearest neighbour distance function G
par(mar=rep(0.2, 4))
plot(cells)
plot(nztrees)
plot(redwood)
par(mar=c(4,4,2,1)+0.2)
plot(Gest(cells))
plot(Gest(nztrees))
plot(Gest(redwood))
## J-function
par(mar=rep(0.2, 4))
plot(cells)
plot(nztrees)
plot(redwood)
par(mar=c(4,4,2,1)+0.2)
plot(Jest(cells))
plot(Jest(nztrees))
plot(Jest(redwood))
par(mfrow=c(1,1),
    mar=c(4,4,2,1)+0.2)

## versions for inhomogeneous patterns
plot(Finhom(japanesepines))
plot(Ginhom(japanesepines))
plot(Jinhom(japanesepines))

## Display F,G,J,K
plot(allstats(swedishpines))

## Multitype patterns
plot(amacrine)
plot(Kcross(amacrine))
plot(Kdot(amacrine))
I <- (marks(amacrine) == "on")
J <- (marks(amacrine) == "off")
plot(Kmulti(amacrine, I, J))

plot(alltypes(amacrine, "K"))

plot(Lcross(amacrine))
plot(Ldot(amacrine))

plot(pcfcross(amacrine))
plot(pcfdot(amacrine))
plot(pcfmulti(amacrine, I, J))

plot(Gcross(amacrine))
plot(Gdot(amacrine))
plot(Gmulti(amacrine, I, J))
plot(alltypes(amacrine, "G"))

plot(Jcross(amacrine))
plot(Jdot(amacrine))
plot(Jmulti(amacrine,I,J))
plot(alltypes(amacrine, "J"))

plot(alltypes(amacrine, "F"))

plot(Iest(amacrine))

plot(markconnect(amacrine))

## Multitype, inhomogeneous
plot(Kcross.inhom(amacrine))
plot(Kdot.inhom(amacrine))
plot(Kmulti.inhom(amacrine, I, J))
plot(Lcross.inhom(amacrine))
plot(Ldot.inhom(amacrine))

plot(pcfcross.inhom(amacrine))
plot(pcfdot.inhom(amacrine))
plot(pcfmulti.inhom(amacrine, I, J))	

## Numerical marks
plot(markcorr(longleaf))
plot(markvario(longleaf))
plot(Emark(longleaf))
plot(Vmark(longleaf))

## Linear networks
plot(chicago)
plot(linearK(chicago))
plot(linearKcross(chicago))
plot(linearKdot(chicago))
plot(linearpcf(chicago))
plot(linearpcfcross(chicago))
plot(linearpcfdot(chicago))

lam <- rep(intensity(unmark(chicago)), npoints(chicago))
A <- split(chicago)$assault
B <- split(chicago)$burglary
lamA <- rep(intensity(A), npoints(A))
lamB <- rep(intensity(B), npoints(B))
plot(linearKinhom(chicago, lam))
plot(linearKcross.inhom(chicago, "assault", "burglary", lamA, lamB))
plot(linearKdot.inhom(chicago, "assault", lamA, lam))
plot(linearpcfinhom(chicago, lam))
plot(linearpcfcross.inhom(chicago, "assault", "burglary", lamA, lamB))
plot(linearpcfdot.inhom(chicago, "assault", lamA, lam))

plot(linearmarkconnect(chicago))
plot(linearmarkequal(chicago))

rm(I,J,fit)

par(opa)
