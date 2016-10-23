#
# anova.mppm.R
#
# $Revision: 1.10 $ $Date: 2016/10/23 10:36:58 $
#

anova.mppm <- local({

  do.gripe <- function(...) warning(paste(...), call.=FALSE)
  dont.gripe <- function(...) NULL
  tests.choices <- c("Chisq", "LRT", "Rao", "score", "F", "Cp")
  tests.avail <- c("Chisq", "LRT", "Rao", "score")
  tests.random  <- c("Chisq", "LRT")
  tests.Gibbs <- c("Chisq", "LRT")
  totalnquad <- function(fit) sum(sapply(quad.mppm(fit), n.quad))
  totalusedquad <- function(fit) with(fit$Fit$moadf, sum(.mpl.SUBSET))
  fmlaString <- function(z) { paste(as.expression(formula(z))) }
##  interString <- function(z) { as.interact(z)$creator }
  
  anova.mppm <- function(object, ..., test=NULL, adjust=TRUE,
                         fine=FALSE, warn=TRUE) {
    gripe <- if(warn) do.gripe else dont.gripe
    argh <- list(...)

    ## trap outmoded usage
    if("override" %in% names(argh)) {
      gripe("Argument 'override' is superseded and was ignored")
      argh <- argh[-which(names(argh) == "override")]
    }
   
    ## list of models
    objex <- append(list(object), argh)

    ## Check each model is an mppm object
    if(!all(sapply(objex, is.mppm)))
      stop(paste("Arguments must all be", sQuote("mppm"), "objects"))

    ## are all models Poisson?
    pois <- all(sapply(objex, is.poisson.mppm))
    gibbs <- !pois

    ## handle anova for a single object
    expandedfrom1 <- FALSE
    if(length(objex) == 1 && gibbs) {
      ## we can't rely on anova.glm in this case
      ## so we have to re-fit explicitly
      Terms <- drop.scope(object)
      if((nT <- length(Terms)) > 0) {
        ## generate models by adding terms sequentially
        objex <- vector(mode="list", length=nT+1)
        for(n in 1L:nT) {
          ## model containing terms 1, ..., n-1
          fmla <- paste(". ~ . - ", paste(Terms[n:nT], collapse=" - "))
          fmla <- as.formula(fmla)
          objex[[n]] <- update(object, fmla)
        }
        ## full model
        objex[[nT+1L]] <- object
        expandedfrom1 <- TRUE
      }
    }

    ## All models fitted using same method?
    Fits <- lapply(objex, getElement, name="Fit")
    fitter <- unique(unlist(lapply(Fits, getElement, name="fitter")))
    if(length(fitter) > 1)
      stop(paste("Models are incompatible;",
                 "they were fitted by different methods (",
                 paste(fitter, collapse=", "), ")" ))

    ## Choice of test
    if(fitter == "glmmPQL") {
      ## anova.lme requires different format of `test' argument
      ## and does not recognise 'dispersion'
      if(is.null(test))
        test <- FALSE
      else {
        test <- match.arg(test, tests.choices)
        if(!(test %in% tests.random))
          stop(paste("Test", dQuote(test),
                     "is not implemented for random effects models"))
        test <- TRUE
      }
    } else if(!is.null(test)) {
      test <- match.arg(test, tests.choices)
      if(!(test %in% tests.avail))
        stop(paste("test=", dQuote(test), "is not yet implemented"),
             call.=FALSE)
      if(!pois && !(test %in% tests.Gibbs))
        stop(paste("test=", dQuote(test),
                   "is only implemented for Poisson models"),
             call.=FALSE)
    }
  

    ## Extract glm fit objects 
    fitz <- lapply(Fits, getElement, name="FIT")

    ## Ensure all models were fitted using GLM, or all were fitted using GAM
    isgam <- sapply(fitz, inherits, what="gam")
    isglm <- sapply(fitz, inherits, what="glm")
    usegam <- any(isgam)
    if(usegam && any(isglm)) {
      gripe("Models were re-fitted with use.gam=TRUE")
      objex <- lapply(objex, update, use.gam=TRUE)
    }

    opt <- list(test=test, dispersion=1)
    if(fitter == "glmmPQL") opt <- list(test=test)

    ## Finally do the appropriate ANOVA
    result <- do.call(anova, append(fitz, opt))
  
    ## Remove approximation-dependent columns if present
    result[, "Resid. Dev"] <- NULL
    ## replace 'residual df' by number of parameters in model
    if("Resid. Df" %in% names(result)) {
      ## count number of quadrature points used in each model
      nq <- totalusedquad(objex[[1L]])
      result[, "Resid. Df"] <- nq - result[, "Resid. Df"]
      names(result)[match("Resid. Df", names(result))] <- "Npar"
    }

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
      for(i in 1L:5L)
        h <- gsub("\n\n", "\n", h)
      if(length(objex) > 1 && length(h) > 1) {
        ## anova(mod1, mod2, ...)
        ## change names of models
        fmlae <- unlist(lapply(objex, fmlaString))
#        intrx <- unlist(lapply(objex, interString))
        h[2L] <- paste("Model",
                      paste0(1L:length(objex), ":"),
                      fmlae,
#                      "\t",
#                      intrx,
                      collapse="\n")
      }
      ## Add explanation if we did the stepwise thing ourselves
      if(expandedfrom1)
        h <- c(h[1L], "Terms added sequentially (first to last)\n", h[-1])
      ## Contract spaces in output if spatstat.options('terse') >= 2
      if(!waxlyrical('space'))
        h <- gsub("\n$", "", h)
      ## Put back
      attr(result, "heading") <- h
    }

    if(adjust && !pois) {
      ## issue warning, if not already given
      if(warn) warn.once("anovaMppmAdjust",
                         "anova.mppm now computes the *adjusted* deviances",
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
              JH <- vcov(bigger, what="all", new.coef=thetaDot, fine=fine)
#              J   <- if(!logi) JH$Sigma else (JH$Sigma1log+JH$Sigma2log)
#              H   <- if(!logi) JH$A1 else JH$Slog
              J <- JH$fisher
              H <- JH$internals$A1
              G   <- H%*%solve(J)%*%H
              if(df == 1) {
                cfac[i] <- H[-injection,-injection]/G[-injection,-injection]
              } else {
                Res <- lapply(subfits(bigger),
                              residuals,
                              type="score",
                              drop=TRUE, 
                              new.coef=thetaDot, dropcoef=TRUE)
                U <- sumcompatible(lapply(Res, integral.msr), names(thetaDot))
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
          result[["Pr(>Chi)"]] <- c(NA, pchisq(abs(AdjDev[-1L]),
                                               df=abs(result$Df[-1L]),
                                               lower.tail=FALSE))
        class(result) <- class(oldresult)
        attr(result, "heading") <- attr(oldresult, "heading")
      }
    }

    return(result)
  }

  sumcompatible <- function(xlist, required) {
    result <- numeric(length(required))
    names(result) <- required
    for(x in xlist) {
      namx <- names(x)
      if(!all(ok <- (namx %in% required)))
        stop(paste("Internal error in sumcompatible:",
                   "list entry", i, "contains unrecognised",
                   ngettext(sum(!ok), "value", "values"),
                   commasep(sQuote(namx[!ok]))),
             call.=FALSE)
      inject <- match(namx, required)
      result[inject] <- result[inject] + x
    }
    return(result)
  }
    
  anova.mppm
})


