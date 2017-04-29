#' temporary file of additional tests
require(spatstat)

#' leverage and influence for logistic fits
fitP <- ppm(swedishpines ~x, method="logi")
fitS <- ppm(swedishpines ~x, Strauss(9), method="logi")
plot(leverage(fitP))
plot(leverage(fitS))
plot(influence(fitP))
plot(influence(fitS))
plot(dfbetas(fitP))
plot(dfbetas(fitS))

#' disconnected linear network
m <- simplenet$m
m[4,5] <- m[5,4] <- m[6,10] <- m[10,6] <- FALSE
L <- linnet(vertices(simplenet), m)
L
summary(L)
Z <- connected(L, what="components")
#' point pattern with no points in one connected component
X <- rpoislpp(lambda=function(x,y) { 10 * (x < 0.5)}, L)
B <- lineardirichlet(X)
plot(B)
summary(B)
