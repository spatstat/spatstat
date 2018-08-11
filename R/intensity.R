#
# intensity.R
#
# Code related to intensity and intensity approximations
#
#  $Revision: 1.20 $ $Date: 2017/06/05 10:31:58 $
#

intensity <- function(X, ...) {
  UseMethod("intensity")
}

intensity.ppp <- function(X, ..., weights=NULL) {
  n <- npoints(X)
  a <- area(Window(X))
  if(is.null(weights)) {
    ## unweighted case - for efficiency
    if(is.multitype(X)) {
      mks <- marks(X)
      answer <- as.vector(table(mks))/a
      names(answer) <- levels(mks)
    } else answer <- n/a
    return(answer)
  }
  ## weighted case 
  if(is.numeric(weights)) {
    check.nvector(weights, n)
  } else if(is.expression(weights)) {
    # evaluate expression in data frame of coordinates and marks
    df <- as.data.frame(X)
    pf <- parent.frame()
    eval.weights <- try(eval(weights, envir=df, enclos=pf))
    if(inherits(eval.weights, "try-error"))
      stop("Unable to evaluate expression for weights", call.=FALSE)
    if(!check.nvector(eval.weights, n, fatal=FALSE, warn=TRUE))
      stop("Result of evaluating the expression for weights has wrong format")
    weights <- eval.weights
  } else stop("Unrecognised format for argument 'weights'")
  ##
  if(is.multitype(X)) {
    mks <- marks(X)
    answer <- as.vector(tapply(weights, mks, sum))/a
    answer[is.na(answer)] <- 0
    names(answer) <- levels(mks)
  } else {
    answer <- sum(weights)/a
  }
  return(answer)
}

intensity.splitppp <- function(X, ..., weights=NULL) {
  if(is.null(weights))
    return(sapply(X, intensity.ppp))
  if(is.expression(weights))
    return(sapply(X, intensity.ppp, weights=weights))
  if(is.numeric(weights)) {
    fsplit <- attr(X, "fsplit")
    n <- length(fsplit)
    check.nvector(weights, n)
    result <- mapply(intensity.ppp, X, weights=split(weights, fsplit))
    result <- simplify2array(result, higher=FALSE)
    return(result)
  }
  stop("Unrecognised format for weights")
}

intensity.ppm <- function(X, ...) {
  if(!identical(valid.ppm(X), TRUE)) {
    warning("Model is invalid - projecting it")
    X <- project.ppm(X)
  }
  if(is.poisson(X)) {
    if(is.stationary(X)) {
      # stationary univariate/multivariate Poisson
      sX <- summary(X, quick="no variances")
      lam <- sX$trend$value
      if(sX$multitype && sX$no.trend) {
        ## trend is ~1; lam should be replicated for each mark
        lev <- levels(marks(data.ppm(X)))
        lam <- rep(lam, length(lev))
        names(lam) <- lev
      }
      return(lam)
    }
    # Nonstationary Poisson
    return(predict(X, ...))
  }
  # Gibbs process
  if(is.multitype(X))
    stop("Not yet implemented for multitype Gibbs processes")
  # Compute first order term
  if(is.stationary(X)) {
    ## activity parameter
    sX <- summary(X, quick="no variances")
    beta <- sX$trend$value
  } else {
    ## activity function (or values of it, depending on '...')
    beta <- predict(X, ...)
  }
  ## apply approximation
  lambda <- PoisSaddle(beta, fitin(X))
  return(lambda)
}

PoisSaddle <- function(beta, fi) {
  ## apply Poisson-Saddlepoint approximation
  ## given first order term and fitted interaction
  stopifnot(inherits(fi, "fii"))
  inte <- as.interact(fi)
  if(identical(inte$family$name, "pairwise"))
    return(PoisSaddlePairwise(beta, fi))
  if(identical(inte$name, "Geyer saturation process"))
    return(PoisSaddleGeyer(beta, fi))
  if(identical(inte$name, "Area-interaction process"))
    return(PoisSaddleArea(beta, fi))
  stop(paste("Intensity approximation is not yet available for",
             inte$name), call.=FALSE)
}

PoisSaddlePairwise <- function(beta, fi) {
  inte <- as.interact(fi)
  Mayer <- inte$Mayer
  if(is.null(Mayer))
    stop(paste("Sorry, not yet implemented for", inte$name))
  # interaction coefficients
  co <- with(fi, coefs[Vnames[!IsOffset]])
  # compute second Mayer cluster integral
  G <- Mayer(co, inte)
  if(is.null(G) || !is.finite(G)) 
    stop("Internal error in computing Mayer cluster integral")
  if(G < 0)
    stop(paste("Unable to apply Poisson-saddlepoint approximation:",
               "Mayer cluster integral is negative"))
  ## solve
  if(is.im(beta)) {
    lambda <- if(G == 0) beta else eval.im(LambertW(G * beta)/G)
  } else {
    lambda <- if(G == 0) beta else (LambertW(G * beta)/G)
    if(length(lambda) == 1) lambda <- unname(lambda)
  }
  return(lambda)
}


# Lambert's W-function

LambertW <- local({

  yexpyminusx <- function(y,x){y*exp(y)-x}

  W <- function(x) {
    result <- rep.int(NA_real_, length(x))
    ok <- is.finite(x) & (x >= 0)
    if(requireNamespace("gsl", quietly=TRUE)) {
      result[ok] <- gsl::lambert_W0(x[ok])
    } else {
      for(i in which(ok))
        result[i] <- uniroot(yexpyminusx, c(0, x[i]), x=x[i])$root
    }
    return(result)
  }

  W
})

PoisSaddleGeyer <- local({

  PoisSaddleGeyer <- function(beta, fi) {
    gamma <- summary(fi)$sensible$param$gamma
    if(gamma == 1) return(beta)
    inte <- as.interact(fi)
    sat <- inte$par$sat
    R   <- inte$par$r
    #' get probability distribution of Geyer statistic under reference model
    z <- Spatstat.Geyer.Nulldist # from sysdata
    if(is.na(m <- match(sat, z$sat)))
      stop(paste("Sorry, the Poisson-saddlepoint approximation",
                 "is not implemented for Geyer models with sat =", sat),
           call.=FALSE)
    probmatrix <- z$prob[[m]]
    maxachievable <- max(which(colSums(probmatrix) > 0)) - 1
    gammarange <- sort(c(1, gamma^maxachievable))
    #' apply approximation
    betavalues <- beta[]
    nvalues <- length(betavalues)
    lambdavalues <- numeric(nvalues)
    for(i in seq_len(nvalues)) {
      beta.i <- betavalues[i]
      ra <- beta.i * gammarange
      lambdavalues[i] <- uniroot(diffapproxGeyer, ra, beta=beta.i,
                                 gamma=gamma, R=R, sat=sat,
                                 probmatrix=probmatrix)$root
    }
    #' return result in same format as 'beta'
    lambda <- beta
    lambda[] <- lambdavalues
    if(length(lambda) == 1) lambda <- unname(lambda)
    return(lambda)
  }

  diffapproxGeyer <- function(lambda, beta, gamma, R, sat, probmatrix) {
    lambda - approxEpoisGeyerT(lambda, beta, gamma, R, sat, probmatrix)
  }
  approxEpoisGeyerT <- function(lambda, beta=1, gamma=1, R=1, sat=1,
                                probmatrix) {
    #' Compute approximation to E_Pois(lambda) Lambda(0,X) for Geyer
    #' ('probmatrix' contains distribution of geyerT(0, Z_n) for each n,
    #' where 'sat' is given, and Z_n is runifdisc(n, radius=2*R).
    possT <- 0:(ncol(probmatrix)-1)
    possN <- 0:(nrow(probmatrix)-1)
    pN <- dpois(possN, lambda * pi * (2*R)^2)
    EgamT <- pN %*% probmatrix %*% (gamma^possT)
    #' assume that, for n > max(possN),
    #' distribution of T is concentrated on T=sat
    EgamT <- EgamT + (gamma^sat) * (1-sum(pN))
    return(beta * EgamT)
  }

  PoisSaddleGeyer
})

PoisSaddleArea <- local({

  PoisSaddleArea <- function(beta, fi) {
    eta <- summary(fi)$sensible$param$eta
    if(eta == 1) return(beta)
    etarange <- range(c(eta^2, 1.1, 0.9))
    inte <- as.interact(fi)
    R   <- inte$par$r
    #' reference distribution of canonical sufficient statistic
    zeroprob <- Spatstat.Area.Zeroprob
    areaquant <- Spatstat.Area.Quantiles
    # expectation of eta^A_n for each n = 0, 1, ....
    EetaAn <- c(1/eta,
                zeroprob + (1-zeroprob) * colMeans((eta^(-areaquant))))
    #' compute approximation
    betavalues <- beta[]
    nvalues <- length(betavalues)
    lambdavalues <- numeric(nvalues)
    for(i in seq_len(nvalues)) {
      beta.i <- betavalues[i]
      ra <- beta.i * etarange
      lambdavalues[i] <- uniroot(diffapproxArea, ra, beta=beta.i,
                                 eta=eta, r=R,
                                 EetaAn=EetaAn)$root
    }
    #' return result in same format as 'beta'
    lambda <- beta
    lambda[] <- lambdavalues
    if(length(lambda) == 1) lambda <- unname(lambda)
    return(lambda)
  }

  diffapproxArea <- function(lambda, beta, eta, r, EetaAn) {
    lambda - approxEpoisArea(lambda, beta, eta, r, EetaAn)
  }

  approxEpoisArea <- function(lambda, beta=1, eta=1, r=1, EetaAn) {
    #' Compute approximation to E_Pois(lambda) Lambda(0,X) for AreaInter
    mu <- lambda * pi * (2*r)^2
    zeta <- pi^2/2 - 1
    theta <-  -log(eta)
    zetatheta <- zeta * theta

    #' contribution from tabulated values
    Nmax <- length(EetaAn) - 1L
    possN <- 0:Nmax
    qN <- dpois(possN, mu)
    # expectation of eta^A when N ~ poisson (truncated)
    EetaA <- sum(qN * EetaAn)

    #' asymptotics for quite large n
    Nbig <- qpois(0.999, mu)
    qn <- 0
    if(Nbig > Nmax) {
      n <- (Nmax+1):Nbig
      #' asymptotic mean uncovered area conditional on this being positive
      mstarn <- (16/((n+3)^2)) * exp(n * (1/4 - log(4/3)))
      ztm <- zetatheta * mstarn
      ok <- (ztm < 1)
      if(!any(ok)) {
        Nbig <- Nmax
        qn <- 0
      } else {
        if(!all(ok)) {
          Nbig <- max(which(!ok)) - 1
          n <- (Nmax+1):Nbig
          ztm <- ztm[1:((Nbig-Nmax)+1)]
        }
        qn <- dpois(n, mu)
        #' asymptotic  probability of complete coverage
        pstarn <- 1 - pmin(1, 3 * (1 + n^2/(16*pi)) * exp(-n/4))
        Estarn <- (1 - ztm)^(-1/zeta)
        EetaA <- EetaA + sum(qn * (pstarn + (1-pstarn) * Estarn))
      }
    }
    #' for very large n, assume complete coverage, so A = 0
    EetaA <- EetaA + 1 - sum(qN) - sum(qn)
    return(beta * eta * EetaA)
  }

  PoisSaddleArea

})
