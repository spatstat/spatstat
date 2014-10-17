#
#   anova.ppm.R
#
#  $Revision: 1.15 $   $Date: 2014/10/17 05:31:31 $
#

anova.ppm <- function(object, ..., test=NULL, adjust=TRUE) {

  if(!is.null(test)) {
    test <- match.arg(test, c("Chisq", "LRT", "Rao", "F", "Cp"))
    if(!(test %in% c("Chisq", "LRT")))
      stop("test=", dQuote(test), "is not yet implemented")
  }

  ## trap outmoded usage
  argh <- list(...)
  if("override" %in% names(argh)) {
    warning("Argument 'override' is superseded and was ignored")
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

  ## all models fitted by MPL or logi?
  fitmethod <- unlist(lapply(objex, getElement, name="method"))
  logi <- all(fitmethod=="logi")
  if(!(all(fitmethod=="mpl")||logi)) 
    stop(paste("Not all models fitted by maximum pseudolikelihood",
               "or logistic regression method;",
               "comparison not possible"))

  ## Extract glmfit objects 
  fitz <- lapply(objex, getglmfit)

  ## Any trivial models? (uniform Poisson)
  trivial <- unlist(lapply(fitz, is.null))
  
  ## force all non-trivial models to be fitted using same method
  ## (all using GLM or all using GAM)
  isgam <- unlist(lapply(fitz, inherits, what="gam"))
  isglm <- unlist(lapply(fitz, inherits, what="glm"))
  usegam <- any(isgam)
  if(usegam && any(isglm)) {
    warning("Some, but not all, models were fitted with use.gam=TRUE;",
            "refitting all models with use.gam=TRUE.")
    objex[isglm] <- lapply(objex[isglm], update.ppm,
                           forcefit=TRUE, use.gam=TRUE)
    fitz[isglm] <- lapply(objex[isglm], getglmfit)   
  }
  
  ## Force any trivial models to be refitted using GLM or GAM
  if(any(trivial)) {
    ## force them to be fitted using glm
    objex[trivial] <- lapply(objex[trivial], update.ppm,
                             forcefit=TRUE, use.gam=usegam)
    fitz[trivial] <- lapply(objex[trivial], getglmfit)
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
    warn.once("anovaAdjust",
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
          warning(paste("Models", i-1, "and", i, "have the same dimension"))
        } else {
          bignames <- names(coef(bigger))
          smallnames <- names(coef(smaller))
          injection <- match(smallnames, bignames)
          if(any(uhoh <- is.na(injection))) {
            warning(paste("Unable to match",
                          ngettext(sum(uhoh), "coefficient", "coefficients"),
                          commasep(sQuote(smallnames[uhoh])),
                          "of model", ismal, 
                          "to coefficients in model", ibig))
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
      warning(paste("Adjustment to composite likelihood does not account for",
                    "irregular trend parameters (covfunargs)",
                    commasep(sQuote(cfa))),
              call.=FALSE)
    }
  }
  return(result)
}
