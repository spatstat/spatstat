if(dev.cur() <= 1) {
  dd <- getOption("device")
  if(is.character(dd)) dd <- get(dd)
  dd()
}

oldpar <- par(ask = interactive() &&
              (.Device %in% c("X11", "GTK", "windows", "Macintosh")))
par(mfrow=c(1,1))
oldoptions <- options(warn = -1)

# 
#######################################################
#

X <- rpoispp(function(x,y) { 1000 * exp(- 4 * x)}, 1000)
plot(X, main="Inhomogeneous Poisson pattern")

fit.hom <- ppm(X ~1, Poisson())
fit.inhom <- ppm(X ~x, Poisson())

diagnose.ppm(fit.inhom, which="marks", type="Pearson",
             main=c("Mark plot",
               "Circles for positive residual mass",
               "Colour for negative residual density"))

par(mfrow=c(1,2))
diagnose.ppm(fit.hom, which="marks", 
             main=c("Wrong model", "(homogeneous Poisson)", "raw residuals"))
diagnose.ppm(fit.inhom, which="marks", 
             main=c("Right model", "(inhomogeneous Poisson)", "raw residuals"))
par(mfrow=c(1,1))

diagnose.ppm(fit.inhom, which="smooth", main="Smoothed residual field")

par(mfrow=c(1,2))
diagnose.ppm(fit.hom, which="smooth",
             main=c("Wrong model", "(homogeneous Poisson)",
                    "Smoothed residual field"))
diagnose.ppm(fit.inhom, which="smooth",
             main=c("Right model", "(inhomogeneous Poisson)",
                    "Smoothed residual field"))

par(mfrow=c(1,1))
diagnose.ppm(fit.inhom, which="x")

par(mfrow=c(1,2))
diagnose.ppm(fit.hom, which="x",
             main=c("Wrong model", "(homogeneous Poisson)",
                    "lurking variable plot for x"))
diagnose.ppm(fit.inhom, which="x",
             main=c("Right model", "(inhomogeneous Poisson)",
                    "lurking variable plot for x"))

par(mfrow=c(1,1))
diagnose.ppm(fit.hom, type="Pearson",main="standard diagnostic plots")

par(mfrow=c(1,2))
diagnose.ppm(fit.hom, main=c("Wrong model", "(homogeneous Poisson)"))
diagnose.ppm(fit.inhom,  main=c("Right model", "(inhomogeneous Poisson)"))
par(mfrow=c(1,1))


# 
#######################################################
#  LEVERAGE/INFLUENCE

plot(leverage(fit.inhom))

plot(influence(fit.inhom))

plot(dfbetas(fit.inhom))

# 
#######################################################
#  COMPENSATORS

## Takes a long time...
CF <- compareFit(listof(hom=fit.hom, inhom=fit.inhom),
                 Kcom, same="iso", different="icom")
plot(CF, main="model compensators", legend=FALSE)
legend("topleft",
       legend=c("empirical K function", "compensator of CSR",
         "compensator of inhomogeneous Poisson"), lty=1:3, col=1:3)

# 
#######################################################
#  Q - Q  PLOTS
#
qqplot.ppm(fit.hom, 40) 
#conclusion: homogeneous Poisson model is not correct
title(main="Q-Q plot of smoothed residuals")

qqplot.ppm(fit.inhom, 40) # TAKES A WHILE...
title(main=c("Right model", "(inhomogeneous Poisson)",
             "Q-Q plot of smoothed residuals"))
# conclusion: fitted inhomogeneous Poisson model looks OK
# 
#######################################################
#
plot(cells)
fitPoisson <- ppm(cells ~1, Poisson())
diagnose.ppm(fitPoisson, 
             main=c("CSR fitted to cells data",
                    "Raw residuals",
                    "No suggestion of departure from CSR"))
diagnose.ppm(fitPoisson, type="pearson",
             main=c("CSR fitted to cells data",
                    "Pearson residuals",
                    "No suggestion of departure from CSR"))
# These diagnostic plots do NOT show evidence of departure from uniform Poisson

plot(Kcom(fitPoisson), cbind(iso, icom) ~ r)
plot(Gcom(fitPoisson), cbind(han, hcom) ~ r)

# K compensator DOES show strong evidence of departure from uniform Poisson

qqplot.ppm(fitPoisson, 40)
title(main=c("CSR fitted to cells data",
        "Q-Q plot of smoothed raw residuals",
        "Strong suggestion of departure from CSR"))
           
# Q-Q plot DOES show strong evidence of departure from uniform Poisson.
#
fitStrauss <- ppm(cells ~1, Strauss(r=0.1))
diagnose.ppm(fitStrauss, 
             main=c("Strauss model fitted to cells data",
                    "Raw residuals"))
diagnose.ppm(fitStrauss, type="pearson",
             main=c("Strauss model fitted to cells data",
                    "Pearson residuals"))

plot(Kcom(fitStrauss), cbind(iso, icom) ~ r)
plot(Gcom(fitStrauss), cbind(han, hcom) ~ r)

# next line takes a LOOONG time ...
qqplot.ppm(fitStrauss, 40, type="pearson")
title(main=c("Strauss model fitted to cells data",
        "Q-Q plot of smoothed Pearson residuals",
        "Suggests adequate fit")) 
# Conclusion: Strauss model seems OK
# 
#######################################################
#
plot(nztrees)
fit <- ppm(nztrees ~1, Poisson())
diagnose.ppm(fit, type="pearson")
title(main=c("CSR fitted to NZ trees",
             "Pearson residuals"))
diagnose.ppm(fit, type="pearson", cumulative=FALSE)
title(main=c("CSR fitted to NZ trees",
             "Pearson residuals (non-cumulative)"))
lurking(fit, expression(x), type="pearson", cumulative=FALSE,
        splineargs=list(spar=0.3))
# Sharp peak at right is suspicious
qqplot.ppm(fit, 40, type="pearson")
title(main=c("CSR fitted to NZ trees",
        "Q-Q plot of smoothed Pearson residuals"))
# Slight suggestion of departure from Poisson at top right of pattern.
par(oldpar)
options(oldoptions)
