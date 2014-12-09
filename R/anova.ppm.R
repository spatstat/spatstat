#
#   anova.ppm.R
#
#  $Revision: 1.17 $   $Date: 2014/12/09 02:39:39 $
#

anova.ppm <- function(object, ..., test=NULL, adjust=TRUE, warn=TRUE) {

  gripe <-
    if(warn) function(...) warning(paste(...), call.=FALSE) else
    function(...) NULL
  
  if(!is.null(test)) {
    test <- match.arg(test, c("Chisq", "LRT", "Rao", "F", "Cp"))
    if(!(test %in% c("Chisq", "LRT")))
      stop("test=", dQuote(test), "is not yet implemented")
  }

  ## trap outmoded usage
  argh <- list(...)
  if("override" %in% names(argh)) {
    gripe("Argument 'override' is superseded and was ignored")
    argh <- argh[-which(names(argh) == "override")]
  }
  
  ## list of models
  objex <- append(list(object), argh)
  if(!all(unlist(lapply(objex, is.ppm))))
    stop(paste("Arguments must all be", sQuote("ppm"), "objects"))

  ## non-Poisson models?
  pois <- all(unlist(lapply(objex, is.poisson.ppm)))

  ## handle anova for a single 'ippm' object  
  expandedfrom1 <- FALSE
  if(length(objex) == 1 && inherits(object, "ippm")) {
    ## we can't rely on anova.glm to get the df right in this case
    ## so we have to re-fit explicitly
    Terms <- drop.scope(object)
    if((nT <- length(Terms)) > 0) {
      ## generate models by adding terms sequentially
      objex <- vector(mode="list", length=nT+1)
      for(n in 1:nT) {
        ## model containing terms 1, ..., n-1
        fmla <- paste(". ~ . - ", paste(Terms[n:nT], collapse=" - "))
        fmla <- as.formula(fmla)
        objex[[n]] <- update(object, fmla)
      }
      ## full model
      objex[[nT+1]] <- object
      expandedfrom1 <- TRUE
    }
  }

  ## all models fitted by same method?
  fitmethod <- unique(unlist(lapply(objex, getElement, name="method")))
  if(length(fitmethod) > 1)
    stop(paste("Models were fitted by different methods",
               commasep(sQuote(fitmethod)), 
               "- comparison is not possible"))
  ## fitted by MPL or logistic?
  if(!(fitmethod %in% c("mpl", "logi")))
    stop(paste("Not implemented for models fitted by method=",
               sQuote(fitmethod)))
  logi <- (fitmethod == "logi")

  refitargs <- list()
  ## fitted to same quadscheme using same edge correction?
  if(length(objex) > 1) {
    ## same data? 
    datas <- lapply(objex, data.ppm)
    samedata <- all(sapply(datas[-1], identical, y=datas[[1]]))
    if(!samedata) stop("Models were fitted to different datasets")
    ## same dummy points?
    quads <- lapply(objex, quad.ppm)
    samequad <- all(sapply(quads[-1], identical, y=quads[[1]]))
    if(!samequad) {
      gripe("Models were re-fitted using a common quadrature scheme")
      sizes <- sapply(quads,
                      function(x) if(inherits(x, "quad")) n.quad(x) else 0)
      imax <- which.max(sizes)
      bigQ <- quads[[imax]]
      refitargs$Q <- bigQ
    }
    ## same edge correction?
    corrxn <- unique(sapply(objex, getElement, name="correction"))
    if(length(corrxn) > 1)
      stop(paste("Models were fitting using different edge corrections",
                 commasep(sQuote(corrxn))))
    if(corrxn == "border") {
      rbord <- unique(sapply(objex, getElement, name="rbord"))
      if(length(rbord) > 1) {
        gripe("Models were re-fitted using a common value of 'rbord'")
        refitargs$rbord <- max(rbord)
      }
    } 
    
    ## Extract glmfit objects 
    fitz <- lapply(objex, getglmfit)

    ## Any trivial models? (uniform Poisson)
    trivial <- unlist(lapply(fitz, is.null))
    if(any(trivial))
      refitargs$forcefit <- TRUE
    
    ## force all non-trivial models to be fitted using same method
    ## (all using GLM or all using GAM)
    isgam <- unlist(lapply(fitz, inherits, what="gam"))
    isglm <- unlist(lapply(fitz, inherits, what="glm"))
    usegam <- any(isgam)
    if(usegam && any(isglm)) {
      gripe("Models were re-fitted with use.gam=TRUE")
      refitargs$use.gam <- TRUE
      refitargs$forcefit <- TRUE
    }

    ## finally refit models
    if(length(refitargs) > 0)
      objex <- do.call(lapply, append(list(X=objex, FUN=update),
                                      refitargs))
  }
  
  ## If any models were fitted by ippm we need to correct the df
  if(any(unlist(lapply(objex, inherits, what="ippm")))) {
    nfree <- unlist(lapply(lapply(objex, logLik), attr, which="df"))
    ncanonical <- unlist(lapply(lapply(objex, coef), length))
    nextra <- nfree - ncanonical
    for(i in seq_along(fitz))
      if(nextra[i] != 0)
        fitz[[i]]$df.residual <- fitz[[i]]$df.residual - nextra[i]
  }
  
  ## Finally do the appropriate ANOVA
  fitz <- lapply(objex, getglmfit)
  result <- do.call("anova", append(fitz, list(test=test, dispersion=1)))

  ## Remove approximation-dependent columns 
  result[, "Resid. Df"] <- NULL
  result[, "Resid. Dev"] <- NULL

  ## edit header 
  if(!is.null(h <- attr(result, "heading"))) {
    ## remove .mpl.Y and .logi.Y from formulae if present
    h <- gsub(".mpl.Y", "", h)
    h <- gsub(".logi.Y", "", h)
    ## delete GLM information if present
    h <- gsub("Model: quasi, link: log", "", h)
    h <- gsub("Model: binomial, link: logit", "", h)
    h <- gsub("Response: ", "", h)
    ## remove blank lines (up to 4 consecutive blanks can occur)
    for(i in 1:5)
      h <- gsub("\n\n", "\n", h)
    if(length(objex) > 1 && length(h) > 1) {
      ## anova(mod1, mod2, ...)
      ## change names of models
      fmlae <- unlist(lapply(objex,
                             function(z) paste(as.expression(formula(z)))))
      intrx <- unlist(lapply(objex, function(z) as.interact(z)$creator))
      h[2] <- paste("Model",
                    paste0(1:length(objex), ":"),
                    fmlae,
                    "\t",
                    intrx,
                    collapse="\n")
    }
    ## Add explanation if we did the stepwise thing ourselves
    if(expandedfrom1)
      h <- c(h[1], "Terms added sequentially (first to last)\n", h[-1])
    ## Contract spaces in output if spatstat.options('terse') >= 2
    if(!waxlyrical('space'))
      h <- gsub("\n$", "", h)
    ## Put back
    attr(result, "heading") <- h
  }
  
  if(adjust && !pois) {
    ## issue warning, if not already given
    if(warn) warn.once("anovaAdjust",
                          "anova.ppm now computes the *adjusted* deviances",
                          "when the models are not Poisson processes.")
    ## Corrected pseudolikelihood ratio 
    nmodels <- length(objex)
    if(nmodels > 1) {
      cfac <- rep(1, nmodels)
      for(i in 2:nmodels) {
        a <- objex[[i-1]]
        b <- objex[[i]]
        df <- length(coef(a)) - length(coef(b))
        if(df > 0) {
          ibig <- i-1
          ismal <- i
        } else {
          ibig <- i
          ismal <- i-1
          df <- -df
        }
        bigger <- objex[[ibig]]
        smaller <- objex[[ismal]]
        if(df == 0) {
          gripe("Models", i-1, "and", i, "have the same dimension")
        } else {
          bignames <- names(coef(bigger))
          smallnames <- names(coef(smaller))
          injection <- match(smallnames, bignames)
          if(any(uhoh <- is.na(injection))) {
            gripe("Unable to match",
                  ngettext(sum(uhoh), "coefficient", "coefficients"),
                  commasep(sQuote(smallnames[uhoh])),
                  "of model", ismal, 
                  "to coefficients in model", ibig)
          } else {
            thetaDot <- 0 * coef(bigger)
            thetaDot[injection] <- coef(smaller)
            JH <- vcov(bigger, what="internals", new.coef=thetaDot)
            J   <- if(!logi) JH$Sigma else (JH$Sigma1log+JH$Sigma2log)
            H   <- if(!logi) JH$A1 else JH$Slog
            G   <- H%*%solve(J)%*%H
            if(df == 1) {
              cfac[i] <- H[-injection,-injection]/G[-injection,-injection]
            } else {
              Res <- residuals(bigger, type="score",
                               new.coef=thetaDot, drop=TRUE)
              U <- integral.msr(Res)
              Uo <- U[-injection]
              Uo <- matrix(Uo, ncol=1)
              Hinv <- solve(H)
              Ginv <- solve(G)
              Hoo <- Hinv[-injection,-injection, drop=FALSE]
              Goo <- Ginv[-injection,-injection, drop=FALSE]
              ScoreStat <- t(Uo) %*% Hoo %*% solve(Goo) %*% Hoo %*% Uo
              cfac[i] <- ScoreStat/(t(Uo) %*% Hoo %*% Uo)
            }
          }
        }
      }
      ## apply Pace et al (2011) adjustment to pseudo-deviances
      ## (save attributes of 'result' for later reinstatement)
      oldresult <- result
      result$Deviance <- AdjDev <- result$Deviance * cfac
      cn <- colnames(result)
      colnames(result)[cn == "Deviance"] <- "AdjDeviance"
      if("Pr(>Chi)" %in% colnames(result)) 
        result[["Pr(>Chi)"]] <- c(NA, pchisq(abs(AdjDev[-1]),
                                       df=abs(result$Df[-1]),
                                       lower.tail=FALSE))
      class(result) <- class(oldresult)
      attr(result, "heading") <- attr(oldresult, "heading")
    }
    if(any(unlist(lapply(objex, inherits, what="ippm")))) {
      ## calculation does not include 'covfunargs'
      cfa <- lapply(lapply(objects, getElement, name="confunargs"), names)
      cfa <- unique(unlist(cfa))
      gripe("Adjustment to composite likelihood does not account for",
            "irregular trend parameters (covfunargs)",
            commasep(sQuote(cfa)))
    }
  }
  return(result)
}
