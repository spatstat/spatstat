#  Variance-covariance matrix for mppm objects
#
# $Revision: 1.15 $ $Date: 2016/02/16 10:08:00 $
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
    collectmom <- (what %in% c("internals", "all"))
    subs <- subfits(object, what="basicmodels")
    n <- length(subs)
    guts <- lapply(subs, vcov, what="internals",
                   matrix.action=matrix.action,
                   gam.action=gam.action,
                   logi.action=logi.action,
                   dropcoef=TRUE,
                   ...)
    fish <- lapply(guts, getElement, name="fisher")
    a1   <- lapply(guts, getElement, name="A1")
    a2   <- lapply(guts, getElement, name="A2")
    a3   <- lapply(guts, getElement, name="A3")
    a1 <- mergeAlternatives(fish, a1)
    cnames <- unique(unlist(lapply(c(a1, a2, a3), colnames)))
    if(collectmom) {
      sufs <- lapply(guts, getElement, name="suff")
      moms <- lapply(guts, getElement, name="mom")
      sufs <- mergeAlternatives(sufs, moms)
      cnames <- unique(c(cnames, unlist(lapply(sufs, colnames))))
    }
    nc <- length(cnames)
    A1 <- A2 <- A3 <- matrix(0, nc, nc, dimnames=list(cnames, cnames))
    if(collectmom)
      Mom <- matrix(, 0, nc, dimnames=list(character(0), cnames))
    for(i in seq_len(n)) {
      coefnames.i <- names(coef(subs[[i]]))
      A1 <- addsubmatrix(A1, a1[[i]], coefnames.i)
      A2 <- addsubmatrix(A2, a2[[i]], coefnames.i)
      A3 <- addsubmatrix(A3, a3[[i]], coefnames.i)
      if(collectmom) Mom <- bindsubmatrix(Mom, sufs[[i]])
    }
    internals <- list(A1=A1, A2=A2, A3=A3)
    if(collectmom)
      internals <- c(internals, list(suff=Mom))
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

  vcov.mppm
  
})
