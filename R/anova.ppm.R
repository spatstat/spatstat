#
#   anova.ppm.R
#
#  $Revision: 1.21 $   $Date: 2015/04/23 14:29:50 $
#

anova.ppm <- local({

  do.gripe <- function(...) warning(paste(...), call.=FALSE)
  dont.gripe <- function(...) NULL
  nquad <- function(x) { if(inherits(x, "quad")) n.quad(x) else 0 }
  fmlaString <- function(z) { paste(as.expression(formula(z))) }
  interString <- function(z) { as.interact(z)$creator }

  anova.ppm <- function(object, ..., test=NULL, adjust=TRUE, warn=TRUE,
                        fine=FALSE) {
    gripe <- if(warn) do.gripe else dont.gripe
    if(!is.null(test)) {
      test <- match.arg(test, c("Chisq", "LRT", "Rao", "score", "F", "Cp"))
      if(test == "score") test <- "Rao"
      if(!(test %in% c("Chisq", "LRT", "Rao")))
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
    if(!all(sapply(objex, is.ppm)))
      stop(paste("Arguments must all be", sQuote("ppm"), "objects"))
    
    ## all models Poisson?
    pois <- all(sapply(objex, is.poisson.ppm))
    gibbs <- !pois
    ## any models fitted by ippm?
    newton <- any(sapply(objex, inherits, what="ippm"))
    
    if(gibbs && !is.null(test) && test == "Rao")
      stop("Score test is only implemented for Poisson models",
           call.=FALSE)
    
    ## handle anova for a single object
    expandedfrom1 <- FALSE
    if(length(objex) == 1 && (gibbs || newton)) {
      ## we can't rely on anova.glm in this case
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
    fitmethod <- unique(sapply(objex, getElement, name="method"))
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
    fitz <- NULL
  
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
        sizes <- sapply(quads, nquad)
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
      trivial <- sapply(fitz, is.null)
      if(any(trivial))
        refitargs$forcefit <- TRUE
    
      ## force all non-trivial models to be fitted using same method
      ## (all using GLM or all using GAM)
      isgam <- sapply(fitz, inherits, what="gam")
      isglm <- sapply(fitz, inherits, what="glm")
      usegam <- any(isgam)
      if(usegam && any(isglm)) {
        gripe("Models were re-fitted with use.gam=TRUE")
        refitargs$use.gam <- TRUE
        refitargs$forcefit <- TRUE
      }

      ## finally refit models
      if(length(refitargs) > 0) {
        objex <- do.call(lapply, append(list(X=objex, FUN=update),
                                        refitargs))
        fitz <- lapply(objex, getglmfit)
      }
    }
  
    ## Ensure GLM/GAM objects all use the same 'subset'
    subz <-  lapply(objex, getglmsubset)
    if(length(unique(subz)) > 1) {
      subsub <- Reduce("&", subz)
      fitz <- lapply(fitz, refittosubset, sub=subsub)
      gripe("Models were re-fitted after discarding quadrature points",
            "that were illegal under some of the models")
    }
  
    ## If any models were fitted by ippm we need to correct the df
    if(newton) {
      nfree <- sapply(lapply(objex, logLik), attr, which="df")
      ncanonical <- sapply(lapply(objex, coef), length)
      nextra <- nfree - ncanonical
      if(is.null(fitz))
        fitz <- lapply(objex, getglmfit)
      for(i in seq_along(fitz))
        if(nextra[i] != 0)
          fitz[[i]]$df.residual <- fitz[[i]]$df.residual - nextra[i]
    }

    ## Finally do the appropriate ANOVA
    if(is.null(fitz)) fitz <- lapply(objex, getglmfit)
    result <- do.call("anova", append(fitz, list(test=test, dispersion=1)))

    ## Remove approximation-dependent columns if present
    result[, "Resid. Dev"] <- NULL
    ## replace 'residual df' by number of parameters in model
    if("Resid. Df" %in% names(result)) {
      ## count number of quadrature points used in each model
      obj1 <- objex[[1]]
      ss <- getglmsubset(obj1)
      nq <- if(!is.null(ss)) sum(ss) else n.quad(quad.ppm(obj1))
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
      for(i in 1:5)
        h <- gsub("\n\n", "\n", h)
      if(length(objex) > 1 && length(h) > 1) {
        ## anova(mod1, mod2, ...)
        ## change names of models
        fmlae <- sapply(objex, fmlaString)
        intrx <- sapply(objex, interString)
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
  
    if(adjust && gibbs) {
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
              JH <- vcov(bigger, what="internals", new.coef=thetaDot, fine=fine)
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
    }

    if(newton) {
      ## calculation does not include 'covfunargs'
      cfa <- lapply(lapply(objex, getElement, name="covfunargs"), names)
      cfa <- unique(unlist(cfa))
      action <- if(adjust && gibbs) "Adjustment to composite likelihood" else
                if(test == "Rao") "Score test calculation" else NULL
      if(!is.null(action)) 
        gripe(action, "does not account for",
              "irregular trend parameters (covfunargs)",
              commasep(sQuote(cfa)))
    }
    return(result)
  }

  refittosubset <- function(fut, sub) {
    etf <- environment(terms(fut))
    gd <- get("glmdata", envir=etf)
    gd$.mpl.SUBSET <- sub
    assign("glmdata", gd, envir=etf)
    up <- update(fut, evaluate=FALSE)
    eval(up, envir=etf)
  }

  anova.ppm
})
