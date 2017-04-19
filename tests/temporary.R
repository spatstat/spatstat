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
