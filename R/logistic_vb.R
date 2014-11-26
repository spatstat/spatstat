#' Variational Bayesian Logistic regression
#' 
#' author: Tuomas Rajala < tuomas.rajala a iki.fi >
#' 
#' Special version for 'spatstat'

####################################################
#' Used inside ppm
vblogit.fmla <- function(formula, offset, data, subset, weights, verbose=FALSE, epsilon=0.01, ...) {
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  offset <- model.offset(mf)
  y <- model.response(mf, "any")
  X <- model.matrix(mt, mf) 
  colnames(X)[1] <- "(Intercept)"
  Vnames <- colnames(X)
  #' then we fit:
  fit <- vblogit(y=y, X=X, offset=offset, verb=verbose, eps=epsilon, ...)
  #'
  names(fit$coefficients) <- names(fit$coef) <- Vnames
  #' add some variables to conform to summary.ppm
  fit$se <- sqrt(diag(as.matrix(fit$S)))
  fit$call <- match.call(expand.dots=FALSE)
  fit$formula <- formula
  fit$method <- "vblogit"
  fit$model <- mf
  fit$terms <- mt
  fit$offset <- offset
  fit$data <- data
  fit$xlevels <- .getXlevels(mt, mf)
  fit
}
###################################################
# the fitting function:
vblogit <- function(y, X, offset, eps=1e-2, m0, S0, S0i, xi0, verb=FALSE, maxiter=1000, ...) {
  ### Logistic regression using JJ96 idea. Ormeron00 notation.
  ## p(y, w, t) = p(y | w) p(w | t) p(t) 
  ##
  ## Y ~ Bern(logit(Xw + offset))
  ## w  ~ N(m0, S0) iid
  ##
  ## "*0" are fixed priors.
  ##
  cat2 <- if(verb) cat else function(...) NULL
  varnames <- colnames(data.frame(as.matrix(X[1:2,])))
  
  ## Write 
  N <- length(y)
  K <- ncol(X)
  #'
  #'
  #' offset
  if(missing('offset')) offset <- 0
  if(length(offset)<N) offset <- rep(offset, N)[1:N]
  #'
  #'
  #' Priors and initial estimates.
  if(missing(S0))S0   <- diag(1e5, K, K)
  if(missing(S0i))S0i <- solve(S0)
  if(missing(m0))m0   <- rep(0, K)
  #' Constants:
  oo2 <- offset^2
  LE_CONST <- as.numeric( -0.5*t(m0)%*%S0i%*%m0 - 0.5*determinant(S0)$mod + sum((y-0.5)*offset) ) 
  Sm0 <- S0i%*%m0
  #' start values for xi:
  if(missing(xi0))xi0   <- rep(4, N) # something positive
  if(length(xi0)!=N) xi0 <- rep(xi0, N)[1:N]
  
  est <- list(m=m0, S=S0, Si=S0i, xi=xi0)
  #'
  #
  ## helper functions needed:
  lambda <- function(x)  -tanh(x/2)/(4*x)
  gamma <- function(x)  x/2 - log(1+exp(x)) + x*tanh(x/2)/4
  ###
  ## loop
  le <- -Inf
  le_hist <- le
  loop <- TRUE
  iter <- 0
  #' initials:
  la <- lambda(xi0)
  L <- diag( x = la )
  Si <- S0i - 2 * t(X)%*%L%*%X
  S <- solve(Si)
  m <- S%*%( t(X)%*%( (y-0.5) + 2*L%*%offset ) + Sm0  )
  #'
  #' Main loop:
  while(loop){
    old <- le
    #' update variational parameters
    xi2 <- diag(  X%*%( S+m%*%t(m) )%*%t(X)+ 2*(X%*%m)%*%t(offset) ) + oo2
    xi <- sqrt(xi2)
    la <- lambda(xi)
    L <- diag( x = la )
    #' update post covariance
    Si <- S0i - 2 * t(X)%*%L%*%X
    S <- solve(Si)
    #' update post mean
    m <- S%*%( t(X)%*%( (y-0.5) + 2*L%*%offset ) + Sm0  )
    #' compute the log evidence
    le <-  as.numeric( 0.5*determinant(S)$mod + sum( gamma(xi) ) + sum(oo2*la) + 0.5*t(m)%*%Si%*%m + LE_CONST)
    #' check convergence 
    devi <- le - old
    if(devi < 0) warning("Log-evidence decreasing, try different starting values for xi.")
    loop <- abs(devi) > eps & (iter<-iter+1) <= maxiter
    le_hist <- c(le_hist, le)
    cat2("diff:", devi, "             \r")
  }
  if(iter == maxiter) warning("Maximum iteration limit reached.")
  cat2("\n")
  ## done. Compile:
  est <- list(m=m, S=S, Si=Si, xi=xi, L=L)
  #' Marginal evidence
  est$logLik <- le
  #' Compute max logLik with the Bernoulli model, this should be what glm gives:
  est$logLik_ML <- as.numeric( t(y)%*%(X%*%m+offset) - sum( log( 1 + exp(X%*%m+offset)) ) )
  #' Max loglik with the approximation
  est$logLik_ML2 <- as.numeric(  t(y)%*%(X%*%m + offset)  + 
                                   t(m)%*%t(X)%*%L%*%X%*%m - 
                                   0.5*sum(X%*%m) + sum(gamma(xi)) +
                                   2*t(offset)%*%L%*%X%*%m + 
                                   t(offset)%*%L%*%offset - 
                                   0.5 * sum(offset)  )
  #' some additional parts, like in glm output
  est$coefficients <- est$m[,1]
  names(est$coefficients) <- varnames
  est$call <- sys.call()
  est$converged <- !(maxiter==iter)
  #' more additional stuff
  est$logp_hist <- le_hist
  est$parameters <- list(eps=eps, maxiter=maxiter)
  est$priors <- list(m=m0, S=S0)
  est$iterations <- iter
  class(est) <- "vblogit"
  ## return
  est
}

###################################################
#' Predict method
predict.vblogit <- function(object, newdata = NULL,
                            type = c("link", "response", "terms"),
                            se.fit = FALSE,
                            dispersion = NULL,
                            terms = NULL,
                            na.action = na.pass, 
                            ...) {
  type <- match.arg(type)
  if(type != "response") stop("type not supported.")
  if(missing(newdata)) {
    stop("not implemented.")
  }
  else{  # newdata
    #' build the new covariate matrix, inspired by predict.lm
    tt <- terms(object)
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.action, 
                     xlev = object$xlevels)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    offset <- rep(0, nrow(X))
    if (!is.null(off.num <- attr(tt, "offset"))) 
      for (i in off.num) offset <- offset + eval(attr(tt, 
                                                      "variables")[[i + 1]], newdata)
    if (!is.null(object$call$offset)) 
      offset <- offset + eval(object$call$offset, newdata)
    #' predict using probit approximation to logit-function
    sigmoid <- function(e) 1/(1+exp(-e))
    mu <- object$m
    S <- object$S
    mua <- as.numeric(X%*%mu)+offset
    s2a <- diag(X%*%S%*%t(X) )
    predictor <- sigmoid( as.numeric( mua/sqrt(1+pi*s2a/8) ) ) 
    names(predictor) <- rownames(X)
  }
  predictor
}


# ###################################################
# print method
print.vblogit <- function(x, ...) {
  cat("Variational Bayes logistic regression fit\n")
  cat("\nCall: ")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\n")
  cat("Log-likelihood:", x$logLik)
  cat("\nConverged:", x$converged)
  cat("\nConvergence threshold:", x$parameters$eps)
  cat("\nIterations / max:", x$iterations, "/", x$parameters$maxiter)
  cat("\n* Caution: the estimates are conditional on convergence.\n")
}
####################################################
# vblogit family method
family.vblogit <- function(x, ...) binomial()

####################################################
#' vblogit fit summary method
summary.vblogit <- function(x, ...) {
  cat("Variational Bayes logistic regression fit\n")
  cat("\nCall: ")
  print(x$call)
  cat("\nCoefficients and posterior 95% central regions:\n")
  vna <- names(x$coefficients)
  s <- sqrt(diag(x$S))
  q0 <- qnorm(c(0.025, 0.975))
  m <- as.numeric(x$m)
  df <- data.frame(estimate=m, "low 0.05"=m+s*q0[1], "high 97.5"=m+s*q0[2], "prior mean"=x$priors$m, "prior var"=diag(x$priors$S))
  rownames(df) <- vna
  print(df)
  cat("\n")
  cat("Lower bound for log-likelihood:", x$logLik)
  cat("\n")
}

####################################################
# Coef
coef.vblogit <- function(x, ...) x$coefficients

####################################################
# Log-evidence
logLik.vblogit <- function(x, ...) {
  x$logLik
}


