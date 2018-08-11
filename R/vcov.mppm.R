#  Variance-covariance matrix for mppm objects
#
# $Revision: 1.19 $ $Date: 2017/12/11 04:10:22 $
#
#

vcov.mppm <- local({

  errhandler <- function(whinge, err) {
    switch(err,
           fatal=stop(whinge),
           warn={
             warning(whinge)
             return(NA)
           },
           null= return(NULL),
           stop(paste("Unrecognised option: err=", dQuote(err))))
  }
    
  vcov.mppm <- function(object, ..., what="vcov", err="fatal") {

    what <- match.arg(what,
                      c("vcov", "corr", "fisher", "Fisher", "internals", "all"))
    if(what == "Fisher") what <- "fisher"

    if(is.poisson.mppm(object) && object$Fit$fitter == "glm") 
      return(vcmPois(object, ..., what=what, err=err))

    return(vcmGibbs(object, ..., what=what, err=err))
  }

  vcmPois <- function(object, ..., what, err) {
    # legacy algorithm for Poisson case
    
    gf <- object$Fit$FIT
    gd <- object$Fit$moadf
    wt <- gd$.mpl.W
    fi <- fitted(gf)

    fo <- object$trend
    if(is.null(fo)) fo <- (~1)

    mof <- model.frame(fo, gd)
    mom <- model.matrix(fo, mof)
    momnames <- dimnames(mom)[[2]]

    fisher <- sumouter(mom, fi * wt)
    dimnames(fisher) <- list(momnames, momnames)

    switch(what,
           fisher = { return(fisher) },
           vcov   = {
             vc <- try(solve(fisher), silent=(err == "null"))
             if(inherits(vc, "try-error"))
               return(errhandler("Fisher information is singular", err))
             else
               return(vc)
           },
           corr={
             co <- try(solve(fisher), silent=(err == "null"))
             if(inherits(co, "try-error"))
               return(errhandler("Fisher information is singular", err))
             sd <- sqrt(diag(co))
             return(co / outer(sd, sd, "*"))
           })
  }

  vcmGibbs <- function(object, ..., what, err,
                       matrix.action=c("warn", "fatal", "silent"),
                       gam.action=c("warn", "fatal", "silent"),
                       logi.action=c("warn", "fatal", "silent")
                       ) {
    if(!missing(err)) {
      if(err == "null") err <- "silent" 
      matrix.action <-
        if(missing(matrix.action)) err else match.arg(matrix.action)
      gam.action <- if(missing(gam.action)) err else match.arg(gam.action)
      logi.action <- if(missing(logi.action)) err else match.arg(logi.action)
    } else {
      matrix.action <- match.arg(matrix.action)
      gam.action <- match.arg(gam.action)
      logi.action <- match.arg(logi.action)
    }
    #' initialise
    cnames <- names(fixed.effects(object))
    nc <- length(cnames)
    A2 <- A3 <- matrix(0, nc, nc, dimnames=list(cnames, cnames))    
    #' (1) Compute matrix A1 directly
    glmdata <- object$Fit$moadf
    glmsub  <- glmdata$.mpl.SUBSET
    wt      <- glmdata$.mpl.W
    mom <- model.matrix(object)
    lam <- unlist(fitted(object))
    A1 <- sumouter(mom, lam * wt * glmsub)
    #' (2) compute A2 and A3 matrices of submodels
    subs <- subfits(object, what="basicmodels")
    n <- length(subs)
    guts <- lapply(subs,
                   vcov,
                   what="internals",
                   matrix.action=matrix.action,
                   gam.action=gam.action,
                   logi.action=logi.action,
                   dropcoef=TRUE,
                   ...)
    a2   <- lapply(guts, getElement, name="A2")
    a3   <- lapply(guts, getElement, name="A3")
    #' (3) map into full model
    #' Identify the (unique) active interaction in each row
    activeinter <- active.interactions(object)
    #' interaction names (in glmdata)
    Vnamelist <- object$Fit$Vnamelist
    #' Each a2[[i]] and a3[[i]] refer to this interaction (eg 'str')
    #' but may contribute to several coefficients of the full model
    #' e.g.  'str' -> str:id -> 'str', 'str:id2'
    #' Determine which canonical variables of full model are active in each row
    mats <- split.data.frame(mom, glmdata$id)
    activevars <- t(sapply(mats, notallzero))
    #' dependence map of canonical variables of full model
    #'     on the original variables/interactions
    md <- model.depends(object$Fit$FIT)
    #' process each row, summing A2 and A3
    for(i in seq_len(n)) {
      #' the submodel in this row
      subi <- subs[[i]]
      #' contributes to second order terms only if non-Poisson
      if(!is.poisson(subi)) {
        cnames.i <- names(coef(subi))
        a2i <- a2[[i]]
        a3i <- a3[[i]]
        #' the (unique) tag name of the interaction in this model
        tagi <- colnames(activeinter)[activeinter[i,]]
        #' the corresponding variable name(s) in glmdata and coef(subi)
        vni <- Vnamelist[[tagi]]
        #' retain only the interaction rows & columns (the rest are zero anyway)
        e <- cnames.i %in% vni
        a2i <- a2i[e, e, drop=FALSE]
        a3i <- a3i[e, e, drop=FALSE]
        cnames.ie <- cnames.i[e]
        #' which coefficients of the full model are active in this row
        acti <- activevars[i,]
        #' for each interaction variable name in the submodel,
        #' find the coefficient(s) in the main model to which it contributes
        nie <- length(cnames.ie)
        cmap <- vector(mode="list", length=nie)
        names(cmap) <- cnames.ie
        for(j in seq_len(nie)) {
          cj <- cnames.ie[j]
          cmap[[j]] <- cnames[ md[,cj] & acti ]
        }
        #' all possible mappings 
        maps <- do.call(expand.grid,
                        append(cmap, list(stringsAsFactors=FALSE)))
        nmaps <- nrow(maps)
        if(nmaps == 0) {
          warning("Internal error: Unable to map submodel to full model")
        } else {
          for(irow in 1:nmaps) {
            for(jcol in 1:nmaps) {
              cmi <- as.character(maps[irow,])
              cmj <- as.character(maps[jcol,])
              if(anyDuplicated(cmi) || anyDuplicated(cmj)) {
                warning("Internal error: duplicated labels in submodel map")
              } else {
                A2[cmi,cmj] <- A2[cmi,cmj] + a2i
                A3[cmi,cmj] <- A3[cmi,cmj] + a2i
              }
            }
          }
        }
      }
    }
    internals <- list(A1=A1, A2=A2, A3=A3)
    if(what %in% c("internals", "all"))
      internals <- c(internals, list(suff=mom))
    if(what %in% c("vcov", "corr", "all")) {
      #' variance-covariance matrix required
      U <- checksolve(A1, matrix.action, , "variance")
      vc <- if(is.null(U)) NULL else (U %*% (A1 + A2 + A3) %*% U)
    }
    out <- switch(what,
                  fisher = A1 + A2 + A3,
                  vcov   = vc,
                  corr   = {
                    if(is.null(vc)) return(NULL)
                    sd <- sqrt(diag(vc))
                    vc / outer(sd, sd, "*")
                  },
                  internals = internals,
                  all = list(internals=internals,
                             fisher=A1+A2+A3,
                             varcov=vc,
                             invgrad=A1)
                  )
    return(out)
  }

  addsubmatrix <- function(A, B, guessnames) {
    if(is.null(B)) return(A)
    if(is.null(colnames(B)) && !missing(guessnames)) {
      if(is.character(guessnames))
        guessnames <- list(guessnames, guessnames)
      if(all(lengths(guessnames) == dim(B)))
        colnames(B) <- guessnames
    }
    if(is.null(colnames(B))) {
      #' unusual
      if(!all(dim(A) == dim(B))) 
        stop("Internal error: no column names, and matrices non-conformable")
      A <- A + B
      return(A)
    }
    j <- match(colnames(B), colnames(A))
    if(anyNA(j)) 
      stop("Internal error: unmatched column name(s)")
    A[j,j] <- A[j,j] + B
    return(A)
  }

  bindsubmatrix <- function(A, B) {
    if(is.null(B)) return(A)
    if(is.null(colnames(B))) {
      if(ncol(A) != ncol(B))
        stop("Internal error: no column names, and matrices non-conformable")
      A <- rbind(A, B)
      return(A)
    }
    j <- match(colnames(B), colnames(A))
    if(anyNA(j))
      stop("Internal error: unmatched column name(s)")
    BB <- matrix(0, nrow(B), ncol(A))
    BB[,j] <- B
    A <- rbind(A, BB)
    return(A)
  }

  mergeAlternatives <- function(A, B) {
    okA <- !sapply(A, is.null)
    okB <- !sapply(B, is.null)
    if(any(override <- !okA & okB))
      A[override] <- B[override]
    return(A)
  }

  notallzero <- function(df) { apply(df != 0, 2, any) }
  
  vcov.mppm
  
})
