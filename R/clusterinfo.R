## lookup table of explicitly-known K functions and pcf
## and algorithms for computing sensible starting parameters

.Spatstat.ClusterModelInfoTable <- 
  list(
       Thomas=list(
         ## Thomas process: old par = (kappa, sigma2) (internally used everywhere)
         ## Thomas process: new par = (kappa, scale) (officially recommended for input/output)
         modelname = "Thomas process", # In modelname field of mincon fv obj.
         descname = "Thomas process", # In desc field of mincon fv obj.
         modelabbrev = "Thomas process", # In fitted obj.
         printmodelname = function(...) "Thomas process", # Used by print.kppm
         parnames = c("kappa", "sigma2"),
         clustargsnames = NULL,
         checkpar = function(par, old = TRUE){
             if(is.null(par))
                 par <- c(kappa=1,scale=1)
             if(any(par<=0))
                 stop("par values must be positive.")
             nam <- check.named.vector(par, c("kappa","sigma2"),
                                       onError="null")
             if(is.null(nam)) {
               check.named.vector(par, c("kappa","scale"))
               names(par)[2L] <- "sigma2"
               par[2L] <- par[2L]^2
             }
             if(!old){
                 names(par)[2L] <- "scale"
                 par[2L] <- sqrt(par[2L])
             }
             return(par)
         },
         checkclustargs = function(margs, old = TRUE) list(),
         resolvedots = function(...){
             ## resolve dots for kppm and friends allowing for old/new par syntax
             dots <- list(...)
             nam <- names(dots)
             out <- list()
             if("ctrl" %in% nam){
                 out$ctrl <- dots$ctrl
             } else{
                 out$ctrl <- dots[nam %in% c("p", "q", "rmin", "rmax")]
             }
             chk <- .Spatstat.ClusterModelInfoTable$Thomas$checkpar
             if(!is.null(dots$startpar)) out$startpar <- chk(dots$startpar)
             return(out)
         },
         # density function for the distance to offspring
         ddist = function(r, scale, ...) {
             2 * pi * r * dnorm(r, 0, scale)/sqrt(2*pi*scale^2)
         },
         ## Practical range of clusters
         range = function(...){
             dots <- list(...)
             par <- dots$par
             # Choose the first of the possible supplied values for scale:
             scale <- c(dots$scale, dots$par[["scale"]], dots$sigma, dots$par[["sigma"]])[1L]
             if(is.null(scale))
                 stop("Argument ", sQuote("scale"), " must be given.")
             thresh <- dots$thresh
             if(!is.null(thresh)){
               ## The squared length of isotropic Gaussian (sigma)
               ## is exponential with mean 2 sigma^2
               rmax <- scale * sqrt(2 * qexp(thresh, lower.tail=FALSE))
               ## old code
               ##  ddist <- .Spatstat.ClusterModelInfoTable$Thomas$ddist
               ##  kernel0 <- clusterkernel("Thomas", scale = scale)(0,0)
               ##  f <- function(r) ddist(r, scale = scale)-thresh*kernel0
               ##  rmax <- uniroot(f, lower = scale, upper = 1000 * scale)$root
             } else{
                 rmax <- 4*scale
             }
             return(rmax)
         },
         kernel = function(par, rvals, ...) {
             scale <- sqrt(par[2L])
             dnorm(rvals, 0, scale)/sqrt(2*pi*scale^2)
         },
         isPCP=TRUE,
         ## K-function
         K = function(par,rvals, ...){
           if(any(par <= 0))
             return(rep.int(Inf, length(rvals)))
           pi*rvals^2+(1-exp(-rvals^2/(4*par[2L])))/par[1L]
         },
         ## pair correlation function
         pcf= function(par,rvals, ...){
           if(any(par <= 0))
             return(rep.int(Inf, length(rvals)))
           1 + exp(-rvals^2/(4 * par[2L]))/(4 * pi * par[1L] * par[2L])
         },
         ## sensible starting parameters
         selfstart = function(X) {
           kappa <- intensity(X)
           sigma2 <- 4 * mean(nndist(X))^2
           c(kappa=kappa, sigma2=sigma2)
         },
         ## meaningful model parameters
         interpret = function(par, lambda) {
           kappa <- par[["kappa"]]
           sigma <- sqrt(par[["sigma2"]])
           mu <- if(is.numeric(lambda) && length(lambda) == 1)
             lambda/kappa else NA
           c(kappa=kappa, sigma=sigma, mu=mu)
         },
         ## Experimental: convert to/from canonical cluster parameters
         tocanonical = function(par) {
           kappa <- par[[1L]]
           sigma2 <- par[[2L]]
           c(strength=1/(kappa * sigma2), scale=sqrt(sigma2))
         },
         tohuman = function(can) {
           strength <- can[[1L]]
           scale <- can[[2L]]
           sigma2 <- scale^2
           c(kappa=1/(strength * sigma2), sigma2=sigma2)
         }
         ),
       ## ...............................................
       MatClust=list(
         ## Matern Cluster process: old par = (kappa, R) (internally used everywhere)
         ## Matern Cluster process: new par = (kappa, scale) (officially recommended for input/output)
         modelname = "Matern cluster process", # In modelname field of mincon fv obj.
         descname = "Matern cluster process", # In desc field of mincon fv obj.
         modelabbrev = "Matern cluster process", # In fitted obj.
         printmodelname = function(...) "Matern cluster process", # Used by print.kppm
         parnames = c("kappa", "R"),
         clustargsnames = NULL,
         checkpar = function(par, old = TRUE){
             if(is.null(par))
                 par <- c(kappa=1,scale=1)
             if(any(par<=0))
                 stop("par values must be positive.")
             nam <- check.named.vector(par, c("kappa","R"), onError="null")
             if(is.null(nam)) {
               check.named.vector(par, c("kappa","scale"))
               names(par)[2L] <- "R"
             }
             if(!old){
                 names(par)[2L] <- "scale"
             }
             return(par)
         },
         # density function for the distance to offspring
         ddist = function(r, scale, ...) {
             ifelse(r>scale, 0, 2 * r / scale^2)
         },
         ## Practical range of clusters
         range = function(...){
             dots <- list(...)
             par <- dots$par
             # Choose the first of the possible supplied values for scale:
             scale <- c(dots$scale, dots$par[["scale"]], dots$R, dots$par[["R"]])[1L]
             if(is.null(scale))
                 stop("Argument ", sQuote("scale"), " must be given.")
           if(!is.null(dots$thresh))
               warning("Argument ", sQuote("thresh"), " is ignored for Matern Cluster model")
             return(scale)
         },
         checkclustargs = function(margs, old = TRUE) list(),
         resolvedots = function(...){
             ## resolve dots for kppm and friends allowing for old/new par syntax
             dots <- list(...)
             nam <- names(dots)
             out <- list()
             if("ctrl" %in% nam){
                 out$ctrl <- dots$ctrl
             } else{
                 out$ctrl <- dots[nam %in% c("p", "q", "rmin", "rmax")]
             }
             chk <- .Spatstat.ClusterModelInfoTable$MatClust$checkpar
             if(!is.null(dots$startpar)) out$startpar <- chk(dots$startpar)
             return(out)
         },
         kernel = function(par, rvals, ...) {
             scale <- par[2L]
             ifelse(rvals>scale, 0, 1/(pi*scale^2))
         },
         isPCP=TRUE,
         K = function(par,rvals, ..., funaux){
           if(any(par <= 0))
             return(rep.int(Inf, length(rvals)))
           kappa <- par[1L]
           R <- par[2L]
           Hfun <- funaux$Hfun
           y <- pi * rvals^2 + (1/kappa) * Hfun(rvals/(2 * R))
           return(y)
         },
         pcf= function(par,rvals, ..., funaux){
             if(any(par <= 0))
               return(rep.int(Inf, length(rvals)))
             kappa <- par[1L]
             R <- par[2L]
             g <- funaux$g
             y <- 1 + (1/(pi * kappa * R^2)) * g(rvals/(2 * R))
             return(y)
           },
         funaux=list(
           Hfun=function(zz) {
             ok <- (zz < 1)
             h <- numeric(length(zz))
             h[!ok] <- 1
             z <- zz[ok]
             h[ok] <- 2 + (1/pi) * (
                                    (8 * z^2 - 4) * acos(z)
                                    - 2 * asin(z)
                                    + 4 * z * sqrt((1 - z^2)^3)
                                    - 6 * z * sqrt(1 - z^2)
                                    )
             return(h)
           },
           DOH=function(zz) {
             ok <- (zz < 1)
             h <- numeric(length(zz))
             h[!ok] <- 0
             z <- zz[ok]
             h[ok] <- (16/pi) * (z * acos(z) - (z^2) * sqrt(1 - z^2))
             return(h)
           },
           ## g(z) = DOH(z)/z has a limit at z=0.
           g=function(zz) {
             ok <- (zz < 1)
             h <- numeric(length(zz))
             h[!ok] <- 0
             z <- zz[ok]
             h[ok] <- (2/pi) * (acos(z) - z * sqrt(1 - z^2))
             return(h)
           }),
         ## sensible starting paramters
         selfstart = function(X) {
           kappa <- intensity(X)
           R <- 2 * mean(nndist(X)) 
           c(kappa=kappa, R=R)
         },
         ## meaningful model parameters
         interpret = function(par, lambda) {
           kappa <- par[["kappa"]]
           R     <- par[["R"]]
           mu    <- if(is.numeric(lambda) && length(lambda) == 1)
             lambda/kappa else NA           
           c(kappa=kappa, R=R, mu=mu)
         }
         ),
       ## ...............................................
       Cauchy=list(
         ## Neyman-Scott with Cauchy clusters: old par = (kappa, eta2) (internally used everywhere)
         ## Neyman-Scott with Cauchy clusters: new par = (kappa, scale) (officially recommended for input/output)
         modelname = "Neyman-Scott process with Cauchy kernel", # In modelname field of mincon fv obj.
         descname = "Neyman-Scott process with Cauchy kernel", # In desc field of mincon fv obj.
         modelabbrev = "Cauchy process", # In fitted obj.
         printmodelname = function(...) "Cauchy process", # Used by print.kppm
         parnames = c("kappa", "eta2"),
         clustargsnames = NULL,
         checkpar = function(par, old = TRUE){
             if(is.null(par))
                 par <- c(kappa=1,scale=1)
             if(any(par<=0))
                 stop("par values must be positive.")
             nam <- check.named.vector(par, c("kappa","eta2"), onError="null")
             if(is.null(nam)) {
                 check.named.vector(par, c("kappa","scale"))
                 names(par)[2L] <- "eta2"
                 par[2L] <- (2*par[2L])^2
             }
             if(!old){
                 names(par)[2L] <- "scale"
                 par[2L] <- sqrt(par[2L])/2
             }
             return(par)
         },
         checkclustargs = function(margs, old = TRUE) list(),
         resolvedots = function(...){
             ## resolve dots for kppm and friends allowing for old/new par syntax
             dots <- list(...)
             nam <- names(dots)
             out <- list()
             if("ctrl" %in% nam){
                 out$ctrl <- dots$ctrl
             } else{
                 out$ctrl <- dots[nam %in% c("p", "q", "rmin", "rmax")]
             }
             chk <- .Spatstat.ClusterModelInfoTable$Cauchy$checkpar
             if(!is.null(dots$startpar)) out$startpar <- chk(dots$startpar)
             return(out)
         },
         # density function for the distance to offspring
         ddist = function(r, scale, ...) {
             r/(scale^2) *  (1 + (r / scale)^2)^(-3/2)
         },
         ## Practical range of clusters
         range = function(...){
             dots <- list(...)
             # Choose the first of the possible supplied values for scale:
             scale <- c(dots$scale, dots$par[["scale"]])[1L]
             if(is.null(scale))
                 stop("Argument ", sQuote("scale"), " must be given.")
             thresh <- dots$thresh %orifnull% 0.01
             ## integral of ddist(r) dr is 1 - (1+(r/scale)^2)^(-1/2)
             ## solve for integral = 1-thresh:
             rmax <- scale * sqrt(1/thresh^2 - 1)
             ## old code
             ## ddist <- .Spatstat.ClusterModelInfoTable$Cauchy$ddist
             ## kernel0 <- clusterkernel("Cauchy", scale = scale)(0,0)
             ## f <- function(r) ddist(r, scale = scale)-thresh*kernel0
             ## rmax <- uniroot(f, lower = scale, upper = 1000 * scale)$root
             return(rmax)
         },
         kernel = function(par, rvals, ...) {
             scale <- sqrt(par[2L])/2
             1/(2*pi*scale^2)*((1 + (rvals/scale)^2)^(-3/2))
         },
         isPCP=TRUE,
         K = function(par,rvals, ...){
           if(any(par <= 0))
             return(rep.int(Inf, length(rvals)))
           pi*rvals^2 + (1 - 1/sqrt(1 + rvals^2/par[2L]))/par[1L]
         },
         pcf= function(par,rvals, ...){
           if(any(par <= 0))
             return(rep.int(Inf, length(rvals)))
           1 + ((1 + rvals^2/par[2L])^(-1.5))/(2 * pi * par[2L] * par[1L])
         },
         selfstart = function(X) {
           kappa <- intensity(X)
           eta2 <- 4 * mean(nndist(X))^2
           c(kappa = kappa, eta2 = eta2)
         },
         ## meaningful model parameters
         interpret = function(par, lambda) {
           kappa <- par[["kappa"]]
           omega <- sqrt(par[["eta2"]])/2
           mu <- if(is.numeric(lambda) && length(lambda) == 1)
             lambda/kappa else NA
           c(kappa=kappa, omega=omega, mu=mu)
         }
         ),
       ## ...............................................
       VarGamma=list(
         ## Neyman-Scott with VarianceGamma/Bessel clusters: old par = (kappa, eta) (internally used everywhere)
         ## Neyman-Scott with VarianceGamma/Bessel clusters: new par = (kappa, scale) (officially recommended for input/output)
         modelname = "Neyman-Scott process with Variance Gamma kernel", # In modelname field of mincon fv obj.
         descname = "Neyman-Scott process with Variance Gamma kernel", # In desc field of mincon fv obj.
         modelabbrev = "Variance Gamma process", # In fitted obj.
         printmodelname = function(obj){ # Used by print.kppm
             paste0("Variance Gamma process (nu=",
                    signif(obj$clustargs[["nu"]], 2), ")")
         },
         parnames = c("kappa", "eta"),
         clustargsnames = "nu",
         checkpar = function(par, old = TRUE){
             if(is.null(par))
                 par <- c(kappa=1,scale=1)
             if(any(par<=0))
                 stop("par values must be positive.")
             nam <- check.named.vector(par, c("kappa","eta"), onError="null")
             if(is.null(nam)) {
               check.named.vector(par, c("kappa","scale"))
               names(par)[2L] <- "eta"
             }
             if(!old) names(par)[2L] <- "scale"
             return(par)
         },
         checkclustargs = function(margs, old = TRUE){
             if(!old)
                 margs <- list(nu=margs$nu.ker)
             return(margs)
         },
         resolvedots = function(...){
             ## resolve dots for kppm and friends allowing for old/new par syntax
             dots <- list(...)
             nam <- names(dots)
             out <- list()
             if("ctrl" %in% nam){
                 out$ctrl <- dots$ctrl
             } else{
                 out$ctrl <- dots[nam %in% c("p", "q", "rmin", "rmax")]
             }
             chk <- .Spatstat.ClusterModelInfoTable$VarGamma$checkpar
             if(!is.null(dots$startpar)) out$startpar <- chk(dots$startpar)
             nu <- dots$nu
             if(is.null(nu)){
                 nu <- try(resolve.vargamma.shape(nu.ker=dots$nu.ker, nu.pcf=dots$nu.pcf)$nu.ker,
                           silent = TRUE)
                 if(inherits(nu, "try-error"))
                     nu <- -1/4
             } else{
                 check.1.real(nu)
                 stopifnot(nu > -1/2)
             }
             out$margs <- list(nu.ker=nu, nu.pcf=2*nu+1)
             out$covmodel <- list(type="Kernel", model="VarGamma", margs=out$margs)
             return(out)
         },
         # density function for the distance to offspring
         ddist = function(r, scale, nu, ...) {
             numer <- ((r/scale)^(nu+1)) * besselK(r/scale, nu)
             numer[r==0] <- 0
             denom <- (2^nu) * scale * gamma(nu + 1)
             numer/denom
         },
         ## Practical range of clusters
         range = function(...){
             dots <- list(...)
             # Choose the first of the possible supplied values for scale:
             scale <- c(dots$scale, dots$par[["scale"]])[1L]
             if(is.null(scale))
                 stop("Argument ", sQuote("scale"), " must be given.")
             # Find value of nu:
             extra <- .Spatstat.ClusterModelInfoTable$VarGamma$resolvedots(...)
             nu <- .Spatstat.ClusterModelInfoTable$VarGamma$checkclustargs(extra$margs, old=FALSE)$nu
             if(is.null(nu))
                 stop("Argument ", sQuote("nu"), " must be given.")
             thresh <- dots$thresh
             if(is.null(thresh))
                 thresh <- .001
             ddist <- .Spatstat.ClusterModelInfoTable$VarGamma$ddist
             f1 <- function(rmx) {
               integrate(ddist, 0, rmx, scale=scale, nu=nu)$value - (1 - thresh)
             }
             f <- Vectorize(f1)
             ## old code
             ## kernel0 <- clusterkernel("VarGamma", scale = scale, nu = nu)(0,0)
             ## f <- function(r) ddist(r, scale = scale, nu = nu) - thresh*kernel0
             rmax <- uniroot(f, lower = scale, upper = 1000 * scale)$root
             return(rmax)
         },
         ## kernel function in polar coordinates (no angular argument).
         kernel = function(par, rvals, ..., margs) {
             scale <- as.numeric(par[2L])
             nu <- margs$nu
             if(is.null(nu))
                 stop("Argument ", sQuote("nu"), " is missing.")
             numer <- ((rvals/scale)^nu) * besselK(rvals/scale, nu)
             numer[rvals==0] <- ifelse(nu>0, 2^(nu-1)*gamma(nu), Inf)
             denom <- pi * (2^(nu+1)) * scale^2 * gamma(nu + 1)
             numer/denom
         },
         isPCP=TRUE,
         K = local({
           ## K function requires integration of pair correlation
           xgx <- function(x, par, nu.pcf) {
             ## x * pcf(x) without check on par values
             numer <- (x/par[2L])^nu.pcf * besselK(x/par[2L], nu.pcf)
             denom <- 2^(nu.pcf+1) * pi * par[2L]^2 * par[1L] * gamma(nu.pcf + 1)
             return(x * (1 + numer/denom))
           }
           vargammaK <- function(par,rvals, ..., margs){
             ## margs = list(.. nu.pcf.. ) 
             if(any(par <= 0))
               return(rep.int(Inf, length(rvals)))
             nu.pcf <- margs$nu.pcf
             out <- numeric(length(rvals))
             ok <- (rvals > 0)
             rvalsok <- rvals[ok]
             outok <- numeric(sum(ok))
             for (i in 1:length(rvalsok))
               outok[i] <- 2 * pi * integrate(xgx,
                                              lower=0, upper=rvalsok[i],
                                              par=par, nu.pcf=nu.pcf)$value
             out[ok] <- outok
             return(out)
           }
           ## Initiated integration in sub-subintervals, but it is unfinished!
           ## vargammaK <- function(par,rvals, ..., margs){
           ##   ## margs = list(.. nu.pcf.. ) 
           ##   if(any(par <= 0))
           ##     return(rep.int(Inf, length(rvals)))
           ##   nu.pcf <- margs$nu.pcf
           ##   out <- numeric(length(rvals))
           ##   out[1L] <- if(rvals[1L] == 0) 0 else 
           ##   integrate(xgx, lower=0, upper=rvals[1L],
           ##             par = par, nu.pcf=nu.pcf)$value
           ##   for (i in 2:length(rvals)) {
           ##     delta <- integrate(xgx,
           ##                        lower=rvals[i-1L], upper=rvals[i],
           ##                        par=par, nu.pcf=nu.pcf)
           ##     out[i]=out[i-1L]+delta$value
           ##   }
           ##   return(out)
           ## }
           vargammaK
           }), ## end of 'local'
         pcf= function(par,rvals, ..., margs){
           ## margs = list(..nu.pcf..)
           if(any(par <= 0))
             return(rep.int(Inf, length(rvals)))
           nu.pcf <- margs$nu.pcf
           sig2 <- 1 / (4 * pi * (par[2L]^2) * nu.pcf * par[1L])
           denom <- 2^(nu.pcf - 1) * gamma(nu.pcf)
           rr <- rvals / par[2L]
           ## Matern correlation function
           fr <- ifelseXB(rr > 0,
                        (rr^nu.pcf) * besselK(rr, nu.pcf) / denom,
                        1)
           return(1 + sig2 * fr)
         },
         parhandler = function(..., nu.ker = -1/4) {
           check.1.real(nu.ker)
           stopifnot(nu.ker > -1/2)
           nu.pcf <- 2 * nu.ker + 1
           return(list(type="Kernel",
                       model="VarGamma",
                       margs=list(nu.ker=nu.ker,
                                  nu.pcf=nu.pcf)))
         },
         ## sensible starting values
         selfstart = function(X) {
           kappa <- intensity(X)
           eta <- 2 * mean(nndist(X))
           c(kappa=kappa, eta=eta)
         },
         ## meaningful model parameters
         interpret = function(par, lambda) {
           kappa <- par[["kappa"]]
           omega <- par[["eta"]]
           mu <- if(is.numeric(lambda) && length(lambda) == 1)
             lambda/kappa else NA
           c(kappa=kappa, omega=omega, mu=mu)
         }
         ),
       ## ...............................................
       LGCP=list(
         ## Log Gaussian Cox process: old par = (sigma2, alpha) (internally used everywhere)
         ## Log Gaussian Cox process: new par = (var, scale) (officially recommended for input/output)
         modelname = "Log-Gaussian Cox process", # In modelname field of mincon fv obj.
         descname = "LGCP", # In desc field of mincon fv obj.
         modelabbrev = "log-Gaussian Cox process", # In fitted obj.
         printmodelname = function(...) "log-Gaussian Cox process", # Used by print.kppm
         parnames = c("sigma2", "alpha"),
         checkpar = function(par, old = TRUE){
             if(is.null(par))
                 par <- c(var=1,scale=1)
             if(any(par<=0))
                 stop("par values must be positive.")
             nam <- check.named.vector(par, c("sigma2","alpha"), onError="null")
             if(is.null(nam)) {
                 check.named.vector(par, c("var","scale"))
                 names(par) <- c("sigma2", "alpha")
             }
             if(!old) names(par) <- c("var", "scale")
             return(par)
         },
         checkclustargs = function(margs, old = TRUE) return(margs),
         resolvedots = function(...){
             ## resolve dots for kppm and friends allowing for old/new par syntax
             dots <- list(...)
             nam <- names(dots)
             out <- list()
             if("ctrl" %in% nam){
                 out$ctrl <- dots$ctrl
             } else{
                 out$ctrl <- dots[nam %in% c("p", "q", "rmin", "rmax")]
             }
             chk <- .Spatstat.ClusterModelInfoTable$LGCP$checkpar
             if(!is.null(dots$startpar)) out$startpar <- chk(dots$startpar)
             cmod <- dots$covmodel
             model <- cmod$model %orifnull% dots$model %orifnull% "exponential"
             margs <- NULL
             if(!identical(model, "exponential")) {
               ## get the 'model generator' 
               modgen <- getRandomFieldsModelGen(model)
               attr(model, "modgen") <- modgen
               if(is.null(cmod)){
                 margsnam <- names(formals(modgen))
                 margsnam <- margsnam[!(margsnam %in% c("var", "scale"))]
                 margs <- dots[nam %in% margsnam]
               } else{
                 margs <- cmod[names(cmod)!="model"]
               }
             }
             if(length(margs)==0)
                 margs <- NULL
             out$margs <- margs
             out$model <- model
             out$covmodel <- list(type="Covariance", model=model, margs=margs)
             return(out)
         },
         isPCP=FALSE,
         ## calls relevant covariance function from RandomFields package
         K = function(par, rvals, ..., model, margs) {
           if(any(par <= 0))
             return(rep.int(Inf, length(rvals)))
           if(model == "exponential") {
             ## For efficiency and to avoid need for RandomFields package
             integrand <- function(r,par,...) 2*pi*r*exp(par[1L]*exp(-r/par[2L]))
           } else {
             kraeverRandomFields()
             integrand <- function(r,par,model,margs) {
               modgen <- attr(model, "modgen")
               if(length(margs) == 0) {
                 mod <- modgen(var=par[1L], scale=par[2L])
               } else {
                 mod <- do.call(modgen,
                                append(list(var=par[1L], scale=par[2L]),
                                       margs))
               }
               2*pi *r *exp(RandomFields::RFcov(model=mod, x=r))
             }
           }
           nr <- length(rvals)
           th <- numeric(nr)
           if(spatstat.options("fastK.lgcp")) {
             ## integrate using Simpson's rule
             fvals <- integrand(r=rvals, par=par, model=model, margs=margs)
             th[1L] <- rvals[1L] * fvals[1L]/2
             if(nr > 1)
               for(i in 2:nr)
                 th[i] <- th[i-1L] +
                   (rvals[i] - rvals[i-1L]) * (fvals[i] + fvals[i-1L])/2
           } else {
             ## integrate using 'integrate'
             th[1L] <- if(rvals[1L] == 0) 0 else 
             integrate(integrand,lower=0,upper=rvals[1L],
                       par=par,model=model,margs=margs)$value
             for (i in 2:length(rvals)) {
               delta <- integrate(integrand,
                                  lower=rvals[i-1L],upper=rvals[i],
                                  par=par,model=model,margs=margs)
               th[i]=th[i-1L]+delta$value
             }
           }
           return(th)
         },
         pcf= function(par, rvals, ..., model, margs) {
           if(any(par <= 0))
             return(rep.int(Inf, length(rvals)))
           if(model == "exponential") {
             ## For efficiency and to avoid need for RandomFields package
             gtheo <- exp(par[1L]*exp(-rvals/par[2L]))
           } else {
             kraeverRandomFields()
             modgen <- attr(model, "modgen")
             if(length(margs) == 0) {
               mod <- modgen(var=par[1L], scale=par[2L])
             } else {
               mod <- do.call(modgen,
                              append(list(var=par[1L], scale=par[2L]),
                                     margs))
             }
             gtheo <- exp(RandomFields::RFcov(model=mod, x=rvals))
           }
           return(gtheo)
         },
         parhandler=function(model = "exponential", ...) {
           if(!is.character(model))
             stop("Covariance function model should be specified by name")
           margs <- c(...)
           if(!identical(model, "exponential")) {
             ## get the 'model generator' 
             modgen <- getRandomFieldsModelGen(model)
             attr(model, "modgen") <- modgen
           }
           return(list(type="Covariance", model=model, margs=margs))
         },
         ## sensible starting values
         selfstart = function(X) {
           alpha <- 2 * mean(nndist(X))
           c(sigma2=1, alpha=alpha)
         },
         ## meaningful model parameters
         interpret = function(par, lambda) {
           sigma2 <- par[["sigma2"]]
           alpha  <- par[["alpha"]]
           mu <- if(is.numeric(lambda) && length(lambda) == 1 && lambda > 0)
             log(lambda) - sigma2/2 else NA
           c(sigma2=sigma2, alpha=alpha, mu=mu)
         }
         )
  )

spatstatClusterModelInfo <- function(name, onlyPCP = FALSE) {
  if(inherits(name, "detpointprocfamily"))
    return(spatstatDPPModelInfo(name))
  if(!is.character(name) || length(name) != 1)
    stop("Argument must be a single character string", call.=FALSE)
  TheTable <- .Spatstat.ClusterModelInfoTable
  nama2 <- names(TheTable)
  if(onlyPCP){
    ok <- sapply(TheTable, getElement, name="isPCP")
    nama2 <- nama2[ok]
  } 
  if(!(name %in% nama2))
    stop(paste(sQuote(name), "is not recognised;",
               "valid names are", commasep(sQuote(nama2))),
         call.=FALSE)
  out <- TheTable[[name]]
  return(out)
}

