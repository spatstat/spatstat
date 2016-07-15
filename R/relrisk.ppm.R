##
##  relrisk.ppm.R
##
##  $Revision: 1.7 $ $Date: 2016/07/15 10:21:26 $
##

relrisk.ppm <- local({

  relrisk.ppm <- function(X, ..., at=c("pixels", "points"),
                          relative=FALSE, se=FALSE, 
                          casecontrol=TRUE, control=1, case,
                          ngrid=NULL, window=NULL) {
    stopifnot(is.ppm(X))
    stopifnot(is.multitype(X))
    control.given <- !missing(control)
    case.given <- !missing(case)
    at <- match.arg(at)
    if(!relative && (control.given || case.given)) {
      aa <- c("control", "case")[c(control.given, case.given)]
      nn <- length(aa)
      warning(paste(ngettext(nn, "Argument", "Arguments"),
                    paste(sQuote(aa), collapse=" and "),
                    ngettext(nn, "was", "were"),
                    "ignored, because relative=FALSE"))
    }
    model <- X
    Y <- data.ppm(model)
    types <- levels(marks(Y))
    ntypes <- length(types)
#    np <- length(coef(model))
    ## compute probabilities or risks
    if(ntypes == 2 && casecontrol) {
      if(control.given || !case.given) {
        stopifnot(length(control) == 1)
        if(is.numeric(control)) {
          icontrol <- control <- as.integer(control)
          stopifnot(control %in% 1:2)
        } else if(is.character(control)) {
          icontrol <- match(control, types)
          if(is.na(icontrol)) stop(paste("No points have mark =", control))
        } else
          stop(paste("Unrecognised format for argument", sQuote("control")))
        if(!case.given)
          icase <- 3 - icontrol
      }
      if(case.given) {
        stopifnot(length(case) == 1)
        if(is.numeric(case)) {
          icase <- case <- as.integer(case)
          stopifnot(case %in% 1:2)
        } else if(is.character(case)) {
          icase <- match(case, types)
          if(is.na(icase)) stop(paste("No points have mark =", case))
        } else stop(paste("Unrecognised format for argument", sQuote("case")))
        if(!control.given) 
          icontrol <- 3 - icase
      }
      switch(at,
             pixels= {
               ## estimate is a single image
               ## compute images of intensities of each mark
               lambda.each <- predict(model, ngrid=ngrid, window=window)
               if(!relative) {
                 ## compute probabilities..
                 ## total intensity (image)
                 lambda.all <- Reduce("+", lambda.each)
                 if(!se) {
                   result <- lambda.each[[icase]]/lambda.all
                   result <- killglitches(result)
                 } else {
                   probs <- lapply(lambda.each, "/", e2=lambda.all)
                   probs <- as.solist(lapply(probs, killglitches))
                   estimate <- probs[[icase]]
                   SE <- SEprobPixels(model, probs)[[icase]]
                   SE <- killglitches(SE)
                   result <- list(estimate=estimate, SE=SE)
                 }
               } else {
                 ## relative risks
                 lambda.ctrl <- lambda.each[[icontrol]]
                 if(!se) {
                   result <- lambda.each[[icase]]/lambda.ctrl
                   result <- killglitches(result)
                 } else {
                   risks <- lapply(lambda.each, "/", e2=lambda.ctrl)
                   risks <- as.solist(lapply(risks, killglitches))
                   estimate <- risks[[icase]]
                   SE <- SErelriskPixels(model, risks, icontrol)[[icase]]
                   SE <- killglitches(SE)
                   result <- list(estimate=estimate, SE=SE)
                 }
               }
             },
             points={
               ## compute intensities of each type
               Ycase <- unmark(Y) %mark% factor(types[icase], levels=types)
               Yctrl <- unmark(Y) %mark% factor(types[icontrol], levels=types)
               lambda.case <- predict(model, locations=Ycase)
               lambda.ctrl <- predict(model, locations=Yctrl)
               if(!relative) {
                 ## compute probabilities
                 ## total intensity
                 lambda.all  <- lambda.case + lambda.ctrl
                 prob.case <- lambda.case/lambda.all
                 if(!se) {
                   result <- prob.case
                 } else {
                   probs <- matrix(, length(prob.case), 2)
                   probs[,icase] <- prob.case
                   probs[,icontrol] <- 1 - prob.case
                   SE <- SEprobPoints(model, probs)[,icase]
                   result <- list(estimate=prob.case, SE=SE)
                 }
               } else {
                 ## compute relative risks
                 risk.case <- lambda.case/lambda.ctrl
                 if(!se) {
                   result <- risk.case
                 } else {
                   risks <- matrix(, length(risk.case), 2)
                   risks[,icase] <- risk.case
                   risks[,icontrol] <- 1
                   SE <- SErelriskPoints(model, risks, icontrol)[,icase]
                   result <- list(estimate=risk.case, SE=SE)
                 }
               }
             })
    } else {
      ## several types
      if(relative) {
        ## need 'control' type
        stopifnot(length(control) == 1)
        if(is.numeric(control)) {
          icontrol <- control <- as.integer(control)
          stopifnot(control %in% 1:ntypes)
        } else if(is.character(control)) {
          icontrol <- match(control, types)
          if(is.na(icontrol)) stop(paste("No points have mark =", control))
        } else
          stop(paste("Unrecognised format for argument", sQuote("control")))
      }
      switch(at,
             pixels={
               ## estimate is a list of images
               ## Compute images of intensities of each type
               lambda.each <- predict(model, ngrid=ngrid, window=window)
               if(!relative) {
                 ## compute probabilities...
                 ## image of total intensity
                 lambda.all <- Reduce("+", lambda.each)
                 probs <- lapply(lambda.each, "/", e2=lambda.all)
                 probs <- as.solist(lapply(probs, killglitches))
                 if(!se) {
                   result <- probs
                 } else {
                   SE <- SEprobPixels(model, probs)
                   SE <- as.solist(lapply(SE, killglitches))
                   result <- list(estimate=probs, SE=SE)
                 }
               } else {
                 ## compute relative risks
                 risks <- lapply(lambda.each, "/",
                                 e2=lambda.each[[icontrol]])
                 risks <- as.solist(lapply(risks, killglitches))
                 if(!se) {
                   result <- risks
                 } else {
                   SE <- SErelriskPixels(model, risks, icontrol)
                   SE <- as.solist(lapply(SE, killglitches))
                   result <- list(estimate=risks, SE=SE)
                 }
               }
             },
             points = {
               ## matrix of intensities of each type at each point
               ## rows=locations, cols=types
               lambda.each <- sapply(types,
                                     predictfortype, 
                                     loc=unmark(Y), model=model, types=types)
               if(!relative) {
                 ## compute probabilities
                 lambda.all <- rowSums(lambda.each)
                 probs <- lambda.each/lambda.all
                 if(!se) {
                   result <- probs
                 } else {
                   SE <- SEprobPoints(model, probs)
                   result <- list(estimate=probs, SE=SE)
                 }
               } else {
                 ## compute relative risks
                 risks <- lambda.each/lambda.each[,icontrol]
                 if(!se) {
                   result <- risks
                 } else {
                   SE <- SErelriskPoints(model, risks, icontrol)
                   result <- list(estimate=risks, SE=SE)
                 }
               }
            })
    }
    return(result)
  }

  modmats <- function(model) {
    # model matrices for data locations for each possible mark
    QM <- quad.ppm(model)
    Y <- QM$data
    QR <- quadscheme.replicated(Y, unmark(Y[FALSE]))
    sourceid <- QR$param$sourceid
    ## canonical covariates 
    mm <- model.matrix(model, Q=QR)
    ## mm is a matrix with one column for canonical covariate
    ## and one row for each marked point in QR.
    mm <- cbind(data.frame(".s"=sourceid, ".m"=marks(QR)), mm)
    ## Split by marks 
    ss <- split(mm, mm$.m)
    ## Reorganise into compatible matrices
    zz <- lapply(ss, reorg)
    return(zz)
  }
  
  reorg <- function(x) {
      z <- x
      rownames(z) <- NULL
      z[x$.s, ] <- z
      return(z[,-(1:2), drop=FALSE])
  }

  SErelriskPoints <- function(model, riskvalues, icontrol) {
    ## riskvalues is a matrix with rows=data locations, cols=types
    types <- colnames(riskvalues)
    ntypes <- length(types)
    ## 
    S.um <- modmats(model)
    S.um <- lapply(S.um, as.matrix)
    ## S.um is a list of matrices, one for each possible type,
    ## each matrix having one row per data location 
    dS.um <- lapply(S.um, "-", e2=S.um[[icontrol]])
    R.um <- mapply("*",
                   dS.um,
                   as.list(as.data.frame(riskvalues)),
                   SIMPLIFY=FALSE)
    ## likewise R.um is a list of matrices
    ##
    vc <- vcov(model)
    VAR <- lapply(R.um, quadform, v=vc)
    VAR <- do.call(cbind, VAR)
    SE <- sqrt(VAR)
    colnames(SE) <- types
    return(SE)
  }

  msubtract <- function(z1, z2) mapply("-", e1=z1, e2=z2, SIMPLIFY=FALSE)

  mmultiply <- function(z1, z2) solapply(z1, "*", e2=z2)
  
  SErelriskPixels <- function(model, riskvalues, icontrol) {
    ## riskvalues is an imlist
    types <- names(riskvalues)
    ntypes <- length(types)
    ## canonical covariates
    S.um <- model.images(model)
    ## S.um is a hyperframe with one column for each mark value
    ## and one row for each canonical covariate
    dS.um <- lapply(S.um, msubtract, 
                    z2=S.um[,icontrol,drop=TRUE])
    R.um <- mapply(mmultiply,
                   z1=dS.um,
                   z2=riskvalues,
                   SIMPLIFY=FALSE)
    VAR <- vector(mode="list", length=ntypes)
    ntypes <- length(types)
    vc <- vcov(model)
    ncoef <- nrow(vc)
    for(type in 1:ntypes) {
      v <- 0
      Rum <- R.um[[type]]
      for(i in 1:ncoef) {
        for(j in 1:ncoef) {
          v <- v + Rum[[i]] * vc[i,j] * Rum[[j]]
        }
      }
      VAR[[type]] <- v
    }
    names(VAR) <- types
    VAR <- as.solist(VAR)
    SE <- as.solist(lapply(VAR, sqrt))
    return(SE)
  }


  SEprobPixels <- function(model, probvalues) {
    ## probvalues is an imlist
    types <- names(probvalues)
    ntypes <- length(types)
    ## canonical covariates
    S.um <- model.images(model)
    ## S.um is a hyperframe with one column for each mark value
    ## and one row for each canonical covariate
    ncoef <- length(coef(model))
    Sbar.u <- vector(mode="list", length=ncoef)
    for(k in 1:ncoef)
      Sbar.u[[k]] <- Reduce("+",
                            mapply("*", e1=S.um[k,,drop=TRUE], e2=probvalues,
                                   SIMPLIFY=FALSE))
    ## Sbar.u is a list of images, one for each canonical covariate
    Sdif.um <- lapply(as.list(S.um), 
                      msubtract,
                      z2=Sbar.u)
    ## Sdif.um is a list of lists of images.
    ##   List of length ntypes,
    ##   each entry being an imlist of length ncoef
    P.um <- mapply(mmultiply,
                   Sdif.um, 
                   probvalues, 
                   SIMPLIFY=FALSE)
    ## P.um is same format as Sdif.um
    vc <- vcov(model)
    ncoef <- nrow(vc)
    VAR <- vector(mode="list", length=ntypes)
    for(m in 1:ntypes) {
      v <- 0
      Pum <- P.um[[m]]
      for(i in 1:ncoef) {
        for(j in 1:ncoef) {
          v <- v + Pum[[i]] * vc[i,j] * Pum[[j]]
        }
      }
      VAR[[m]] <- v
    }
    names(VAR) <- types
    VAR <- as.solist(VAR)
    SE <- as.solist(lapply(VAR, sqrt))
  }
  
  SEprobPoints <- function(model, probvalues) {
    ## probvalues is a matrix with row=location and column=type
    types <- colnames(probvalues)
    ntypes <- length(types)
    ## canonical covariates
    S.um <- modmats(model)
    S.um <- lapply(S.um, as.matrix)
    ## S.um is a list of matrices, one for each possible type,
    ## each matrix having rows=locations and cols=covariates
    ## Weight each matrix by its mark probabilities
    SW <- mapply("*",
                 e1=S.um,
                 e2=as.list(as.data.frame(probvalues)),
                 SIMPLIFY=FALSE)
    ## average them
    Sbar.u <- Reduce("+", SW)
    ## Sbar.u is a matrix with rows=locations and cols=covariates
    Sdif.um <- lapply(S.um, "-", e2=Sbar.u)
    ## Sdif.um is a list of matrices like S.um
    P.um <- mapply("*",
                   e1=Sdif.um, 
                   e2=as.list(as.data.frame(probvalues)),
                   SIMPLIFY=FALSE)
    ## P.um likewise
    vc <- vcov(model)
    VAR <- lapply(P.um, quadform, v=vc)
    VAR <- do.call(cbind, VAR)
    SE <- sqrt(VAR)
    colnames(SE) <- types
    return(SE)
  }
  
  predictfortype <- function(type, model, types, loc) {
    predict(model, locations=loc %mark% factor(type, levels=types))
  }

  killglitches <- function(z, eps=.Machine$double.eps) {
    ra <- range(z, finite=TRUE)
    if(max(abs(ra)) < eps) {
      z[] <- 0
      return(z)
    }
    if(diff(ra) < eps) 
      z[] <- mean(z, na.rm=TRUE)
    return(z)
  }

  relrisk.ppm
})

