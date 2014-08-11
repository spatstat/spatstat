#
#
#   rmhmodel.R
#
#   $Revision: 1.59 $  $Date: 2013/04/25 06:37:43 $
#
#

rmhmodel <- function(...) {
  UseMethod("rmhmodel")
}


rmhmodel.rmhmodel <- function(model, ...) {
  # Check for outdated internal format
  # C.par was replaced by C.beta and C.ipar in spatstat 1.22-3 
  if(outdated <- !is.null(model$C.par))
    warning("Outdated internal format of rmhmodel object; rebuilding it")
  if(outdated || (length(list(...)) > 0))
    model <- rmhmodel.list(unclass(model), ...)
  return(model)
}

rmhmodel.list <- function(model, ...) {
  argnames <- c("cif","par","w","trend","types")
  ok <- argnames %in% names(model)
  do.call("rmhmodel.default",
          resolve.defaults(list(...), model[argnames[ok]]))
}

rmhmodel.default <- function(...,
                             cif=NULL, par=NULL, w=NULL, trend=NULL, types=NULL)
{
  extractsecret <- function(..., stopinvalid=TRUE) {
    if(length(list(...)) > 0)
      stop(paste("rmhmodel.default: syntax should be",
                 "rmhmodel(cif, par, w, trend, types)",
                 "with arguments given by name if they are present"),
           call.=FALSE)
    return(list(stopinvalid=stopinvalid))
  }
  stopinvalid <- extractsecret(...)$stopinvalid
  
  # Validate parameters
  if(is.null(cif)) stop("cif is missing or NULL")
  if(is.null(par)) stop("par is missing or NULL")

  if(!is.null(w))
    w <- as.owin(w)
  
  if(!is.character(cif))
    stop("cif should be a character string")

  betamultiplier <- 1

  Ncif <- length(cif)
  if(Ncif > 1) {
    # hybrid
    # check for Poisson components
    ispois <- (cif == 'poisson')
    if(any(ispois)) {
      # validate Poisson components
      Npois <- sum(ispois)
      poismodels <- vector(mode="list", length=Npois)
      parpois <- par[ispois]
      for(i in 1:Npois)
        poismodels[[i]] <- rmhmodel(cif='poisson', par=parpois[[i]],
                                    w=w, trend=NULL, types=types,
                                    stopinvalid=FALSE)
      # consolidate Poisson intensity parameters
      poisbetalist <- lapply(poismodels, function(x){x$C.beta})
      poisbeta <- Reduce("*", poisbetalist)
      if(all(ispois)) {
        # model collapses to a Poisson process
        cif <- 'poisson'
        Ncif <- 1
        par <- list(beta=poisbeta)
        betamultiplier <- 1
      } else {
        # remove Poisson components
        cif <- cif[!ispois]
        Ncif <- sum(!ispois)
        par <- par[!ispois]
        if(Ncif == 1) # revert to single-cif format
          par <- par[[1]]
        # absorb beta parameters 
        betamultiplier <- poisbeta
      }
    }
  }
  
  if(Ncif > 1) {
    # genuine hybrid 
    models <- vector(mode="list", length=Ncif)
    check <- vector(mode="list", length=Ncif)
    for(i in 1:Ncif) 
      models[[i]] <- rmhmodel(cif=cif[i], par=par[[i]],
                              w=w, trend=NULL, types=types,
                              stopinvalid=FALSE)
    C.id  <- unlist(lapply(models, function(x){x$C.id}))
    C.betalist <- lapply(models, function(x){x$C.beta})
    C.iparlist <- lapply(models, function(x){x$C.ipar})
    # absorb beta multiplier into beta parameter of first component
    C.betalist[[1]] <- C.betalist[[1]] * betamultiplier
    # concatenate for use in C
    C.beta     <- unlist(C.betalist)
    C.ipar     <- unlist(C.iparlist)
    check <- lapply(models, function(x){x$check})
    maxr <- max(unlist(lapply(models, function(x){x$reach})))
    ismulti <- unlist(lapply(models, function(x){x$multitype.interact}))
    multi <- any(ismulti)
    # determine whether model exists
    integ <- unlist(lapply(models, function(x) { x$integrable }))
    stabi <- unlist(lapply(models, function(x) { x$stabilising }))
    integrable <- all(integ) || any(stabi)
    stabilising <- any(stabi)
    # string explanations of conditions for validity
    integ.ex <- unlist(lapply(models, function(x){ x$explainvalid$integrable }))
    stabi.ex <- unlist(lapply(models, function(x){ x$explainvalid$stabilising}))
    stabi.oper <- !(stabi.ex %in% c("TRUE", "FALSE"))
    integ.oper <- !(integ.ex %in% c("TRUE", "FALSE"))
    compnames <- if(!any(duplicated(C.id))) paste("cif", sQuote(C.id)) else 
         paste("component", 1:Ncif, paren(sQuote(C.id)))
    if(!integrable && stopinvalid) {
      # model is not integrable: explain why
      ifail <- !integ & integ.oper
      ireason <- paste(compnames[ifail], "should satisfy",
                       paren(integ.ex[ifail], "{"))
      ireason <- verbalogic(ireason, "and")
      if(sum(ifail) <= 1) {
        # There's only one offending cif, so stability is redundant
        sreason <- "FALSE"
      } else {
        sfail <- !stabi & stabi.oper
        sreason <- paste(compnames[sfail], "should satisfy",
                         paren(stabi.ex[sfail], "{"))
        sreason <- verbalogic(sreason, "or")
      }
      reason <- verbalogic(c(ireason, sreason), "or")
      stop(paste("rmhmodel: hybrid model is not integrable; ", reason),
           call.=FALSE)
    } else {
      # construct strings summarising conditions for validity
      if(!any(integ.oper))
        ireason <- as.character(integrable)
      else {
        ireason <- paste(compnames[integ.oper], "should satisfy",
                         paren(integ.ex[integ.oper], "{"))
        ireason <- verbalogic(ireason, "and")
      }
      if(!any(stabi.oper))
        sreason <- as.character(stabilising)
      else {
        sreason <- paste(compnames[stabi.oper], "should satisfy",
                         paren(stabi.ex[stabi.oper], "{"))
        sreason <- verbalogic(sreason, "or")
      }
      ireason <- verbalogic(c(ireason, sreason), "or")
      explainvalid <- list(integrable=ireason,
                           stabilising=sreason)
    }

    out <- list(cif=cif,
                par=par,
                w=w,
                trend=trend,
                types=types,
                C.id=C.id,
                C.beta=C.beta,
                C.ipar=C.ipar,
                C.betalist=C.betalist,
                C.iparlist=C.iparlist,
                check=check,
                multitype.interact=multi,
                integrable=integrable,
                stabilising=stabilising,
                explainvalid=explainvalid,
                reach=maxr)
    class(out) <- c("rmhmodel", class(out))
    return(out)
  }

  # non-hybrid
  
  # Check that this is a recognised model
  # and look up the rules for this model
  rules <- spatstatRmhInfo(cif)
  
  # Map the name of the cif from R to C
  #      (the names are normally identical in R and C,
  #      except "poisson" -> NA)
  C.id <- rules$C.id
  
  # Check that the C name is recognised in C 
  if(!is.na(C.id)) {
    z <- .C("knownCif",
            cifname=as.character(C.id),
            answer=as.integer(0),
            PACKAGE="spatstat")
    ok <- as.logical(z$answer)
    if(!ok)
      stop(paste("Internal error: the cif", sQuote(C.id),
                 "is not recognised in the C code"))
  }

  # Validate the model parameters and reformat them 
  check <- rules$parhandler
  checkedpar <-
    if(!rules$multitype)
      check(par)
    else if(!is.null(types))
      check(par, types)
    else 
      # types vector not given - defer checking
      NULL

  if(!is.null(checkedpar)) {
    stopifnot(is.list(checkedpar))
    stopifnot(!is.null(names(checkedpar)) && all(nzchar(names(checkedpar))))
    stopifnot(names(checkedpar)[[1]] == "beta")
    C.beta  <- unlist(checkedpar[[1]])
    C.beta <- C.beta * betamultiplier
    C.ipar <- as.numeric(unlist(checkedpar[-1]))
  } else {
    C.beta <- C.ipar <- NULL
  }
  
  # Determine whether model is integrable
  integrable <- rules$validity(par, "integrable")
  explainvalid  <- rules$explainvalid

  if(!integrable && stopinvalid) 
    stop(paste("rmhmodel: the model is not integrable; it should satisfy",
               explainvalid$integrable),
         call.=FALSE)
  
  # Determine whether cif is stabilising
  # (i.e. any hybrid including this cif will be integrable)
  stabilising <- rules$validity(par, "stabilising")

  # Calculate reach of model
  mreach <- rules$reach(par)

  
###################################################################
# return augmented list  
  out <- list(cif=cif,
              par=par,
              w=w,
              trend=trend,
              types=types,
              C.id=C.id,
              C.beta=C.beta,
              C.ipar=C.ipar,
              check= if(is.null(C.ipar)) check else NULL,
              multitype.interact=rules$multitype,
              integrable=integrable,
              stabilising=stabilising,
              explainvalid=explainvalid,
              reach=mreach
              )
  class(out) <- c("rmhmodel", class(out))
  return(out)
}

print.rmhmodel <- function(x, ...) {
  verifyclass(x, "rmhmodel")

  cat("Metropolis-Hastings algorithm, model parameters\n")

  Ncif <- length(x$cif)
  cat(paste("Conditional intensity:",
            if(Ncif == 1) "cif=" else "hybrid of cifs",
            commasep(sQuote(x$cif)), "\n"))
  
  if(!is.null(x$types)) {
    if(length(x$types) == 1)
      cat("Univariate process.\n")
    else {
      cat("Multitype process with types =\n")
      print(x$types)
      if(!x$multitype.interact)
        cat("Interaction does not depend on type\n")
    }
  } else if(x$multitype.interact) 
    cat("Multitype process, types not yet specified.\n")
  
  cat("Numerical parameters: par =\n")
  print(x$par)
  if(is.null(x$C.ipar))
    cat("Parameters have not yet been checked for compatibility with types.\n")
  if(is.owin(x$w)) print(x$w) else cat("Window: not specified.\n")
  cat("Trend: ")
  if(!is.null(x$trend)) print(x$trend) else cat("none.\n")
  if(!is.null(x$integrable) && !x$integrable) {
    cat("\n*Warning: model is not integrable and cannot be simulated*\n")
  }
  invisible(NULL)
}

reach.rmhmodel <- function(x, ...) {
  if(length(list(...)) == 0)
    return(x$reach)
  # reach must be recomputed 
  cif <- x$cif
  Ncif <- length(cif)
  pars <- if(Ncif == 1) list(x$par) else x$par
  maxr <- 0
  for(i in seq_len(Ncif)) {
    cif.i <- cif[i]
    par.i <- pars[[i]]
    rules <- spatstatRmhInfo(cif.i)
    rchfun  <- rules$reach
    if(!is.function(rchfun))
      stop(paste("Internal error: reach is unknown for cif=", sQuote(cif.i)),
           call.=FALSE)
    r.i <- rchfun(par.i, ...)
    maxr <- max(maxr, r.i, na.rm=TRUE)
  }
  return(maxr)
}

is.poisson.rmhmodel <- function(x) {
  verifyclass(x, "rmhmodel")
  identical(x$cif, 'poisson')
}

is.stationary.rmhmodel <- function(x) {
  verifyclass(x, "rmhmodel")
  tren <- x$trend
  return(is.null(tren) || is.numeric(tren))
}

as.owin.rmhmodel <- function(W, ..., fatal=FALSE) {
  # W is the rmhmodel object. It contains a window w
  ans <- W$w
  if(is.owin(ans)) return(ans)
  if(fatal) stop("rmhmodel object does not contain a window")
  return(NULL)
}

is.expandable.rmhmodel <- function(x) {
  tren <- x$tren
  ok <- function(z) { is.null(z) || is.numeric(z) || is.function(z) }
  return(if(!is.list(tren)) ok(tren) else all(unlist(lapply(tren, ok))))
}

  
#####  Table of rules for handling rmh models ##################

spatstatRmhInfo <- function(cifname) {
  rules <- .Spatstat.RmhTable[[cifname]]
  if(is.null(rules))
    stop(paste("Unrecognised cif:", sQuote(cifname)), call.=FALSE)
  return(rules)
}
  
.Spatstat.RmhTable <-
  list(
#
# 0. Poisson (special case)
#
       'poisson'=
       list(
            C.id=NA,
            multitype=FALSE,
            parhandler=function(par, ...) {
              ctxt <- "For the Poisson process"
              with(par, forbidNA(beta, ctxt))
              par <- check.named.list(par, "beta", ctxt)
              with(par, explain.ifnot(all(beta >= 0), ctxt))
              return(par)
            },
            validity=function(par, kind) {
              switch(kind,
                     integrable=TRUE,
                     stabilising=FALSE)
            },
            explainvalid=list(integrable="TRUE",stabilising="FALSE"),
            reach = function(par, ...) { return(0) }
            ),
#       
# 1. Strauss.
#       
       'strauss'=
       list(
            C.id="strauss",
            multitype=FALSE,
            parhandler=function(par, ...) {
              ctxt <- "For the strauss cif"
              par <- check.named.list(par, c("beta","gamma","r"), ctxt)
              # treat r=NA as absence of interaction
              par <- within(par, if(is.na(r)) { r <- 0; gamma <- 1 })
              with(par, check.finite(beta, ctxt))
              with(par, check.finite(gamma, ctxt))
              with(par, check.finite(r, ctxt))
              with(par, check.1.real(gamma, ctxt))
              with(par, check.1.real(r,     ctxt))
              with(par, explain.ifnot(all(beta >= 0), ctxt))
              with(par, explain.ifnot(gamma >= 0, ctxt))
              with(par, explain.ifnot(r >= 0, ctxt))
              return(par)
            },
            validity=function(par, kind) {
              gamma <- par$gamma
              switch(kind,
                     integrable=(gamma <= 1),
                     stabilising=(gamma == 0)
                     )
            },
            explainvalid=list(
              integrable="gamma <= 1",
              stabilising="gamma == 0"),
            reach = function(par, ...) {
              r <- par[["r"]]
              g <- par[["gamma"]]
              return(if(g == 1) 0 else r)
            }
            ),
#       
# 2. Strauss with hardcore.
#       
       'straush' =
       list(
            C.id="straush",
            multitype=FALSE,
            parhandler=function(par, ...) {
              ctxt <- "For the straush cif"
              par <- check.named.list(par, c("beta","gamma","r","hc"), ctxt)
              # treat hc=NA as absence of hard core
              par <- within(par, if(is.na(hc)) { hc <- 0 } )
              # treat r=NA as absence of interaction
              par <- within(par, if(is.na(r)) { r <- hc; gamma <- 1 } )
              with(par, check.finite(beta, ctxt))
              with(par, check.finite(gamma, ctxt))
              with(par, check.finite(r, ctxt))
              with(par, check.finite(hc, ctxt))
              with(par, check.1.real(gamma, ctxt))
              with(par, check.1.real(r,     ctxt))
              with(par, check.1.real(hc,     ctxt))
              with(par, explain.ifnot(all(beta >= 0), ctxt))
              with(par, explain.ifnot(gamma >= 0, ctxt))
              with(par, explain.ifnot(r >= 0, ctxt))
              with(par, explain.ifnot(hc >= 0, ctxt))
              with(par, explain.ifnot(hc <= r, ctxt))
              return(par)
            },
            validity=function(par, kind) {
              hc <- par$hc
              gamma <- par$gamma
              switch(kind,
                     integrable=(hc > 0 || gamma <= 1),
                     stabilising=(hc > 0)
                   )
            },
            explainvalid=list(
              integrable="hc > 0 or gamma <= 1",
              stabilising="hc > 0"),
            reach = function(par, ...) {
              h <- par[["hc"]]
              r <- par[["r"]]
              g <- par[["gamma"]]
              return(if(g == 1) h else r)
            }
            ),
#       
# 3. Softcore.
#
       'sftcr' =
       list(
            C.id="sftcr",
            multitype=FALSE,
            parhandler=function(par, ...) {
              ctxt <- "For the sftcr cif"
              par <- check.named.list(par, c("beta","sigma","kappa"), ctxt)
              with(par, check.finite(beta, ctxt))
              with(par, check.finite(sigma, ctxt))
              with(par, check.finite(kappa, ctxt))
              with(par, check.1.real(sigma, ctxt))
              with(par, check.1.real(kappa, ctxt))
              with(par, explain.ifnot(all(beta >= 0), ctxt))
              with(par, explain.ifnot(sigma >= 0, ctxt))
              with(par, explain.ifnot(kappa >= 0 && kappa <= 1, ctxt))
              return(par)
            },
            validity=function(par, kind) {
              switch(kind,
                     integrable=TRUE,
                     stabilising=FALSE)
            },
            explainvalid=list(integrable="TRUE",stabilising="FALSE"),
            reach = function(par, ..., epsilon=0) {
              if(epsilon==0)
                return(Inf)
              kappa <- par[["kappa"]]
              sigma <- par[["sigma"]]
              return(sigma/(epsilon^(kappa/2)))
            }
            ),
#       
# 4. Multitype Strauss.
#       
       'straussm' =
       list(
            C.id="straussm",
            multitype=TRUE,
            parhandler=function(par, types) {
              ctxt <- "For the straussm cif"
              par <- check.named.list(par, c("beta","gamma","radii"), ctxt)

              beta <- par$beta
              gamma <- par$gamma
              r <- par$radii
              ntypes <- length(types)

              check.finite(beta, ctxt)
              check.nvector(beta, ntypes, TRUE, "types")

              MultiPair.checkmatrix(gamma, ntypes, "par$gamma")
	      gamma[is.na(gamma)] <- 1
              check.finite(gamma, ctxt)

              MultiPair.checkmatrix(r, ntypes, "par$radii")
              if(any(nar <- is.na(r))) {
                r[nar] <- 0
                gamma[nar] <- 1
              }
              check.finite(r, ctxt)

              explain.ifnot(all(beta >= 0), ctxt)
              explain.ifnot(all(gamma >= 0), ctxt)
              explain.ifnot(all(r >= 0), ctxt)

              par <- list(beta=beta, gamma=gamma, r=r)
              return(par)
            }, 
            validity=function(par, kind) {
              gamma <- par$gamma
              radii <- par$radii
              dg <- diag(gamma)
              dr <- diag(radii)
              hard <-!is.na(dg) & (dg == 0) & !is.na(dr) & (dr > 0)
              operative <- !is.na(gamma) & !is.na(radii) & (radii > 0)
              switch(kind,
                     stabilising=all(hard),
                     integrable=all(hard) || all(gamma[operative] <= 1))
            },
            explainvalid=list(
              integrable=paste(
                "gamma[i,j] <= 1 for all i and j,",
                "or gamma[i,i] = 0 for all i"),
              stabilising="gamma[i,i] = 0 for all i"),
            reach = function(par, ...) {
              r <- par$radii
              g <- par$gamma
              operative <- ! (is.na(r) | (g == 1))
              return(max(0, r[operative]))
            }
            ),
#       
# 5. Multitype Strauss with hardcore.
#       
       'straushm' = 
       list(
            C.id="straushm",
            multitype=TRUE,
            parhandler=function(par, types) {
              ctxt="For the straushm cif"
              par <- check.named.list(par,
                                      c("beta","gamma","iradii","hradii"),
                                      ctxt)

              beta <- par$beta
              gamma <- par$gamma
              iradii <- par$iradii
              hradii <- par$hradii
              ntypes <- length(types)

              check.nvector(beta, ntypes, TRUE, "types")
              check.finite(beta, ctxt)
              
              MultiPair.checkmatrix(gamma, ntypes, "par$gamma")
              gamma[is.na(gamma)] <- 1
              check.finite(gamma, ctxt)

              MultiPair.checkmatrix(iradii, ntypes, "par$iradii")
              if(any(nar <- is.na(iradii))) {
                iradii[nar] <- 0
                gamma[nar] <- 1
              }
              check.finite(iradii, ctxt)

              MultiPair.checkmatrix(hradii, ntypes, "par$hradii")
              hradii[is.na(hradii)] <- 0
              check.finite(hradii, ctxt)

              explain.ifnot(all(beta >= 0), ctxt)
              explain.ifnot(all(gamma >= 0), ctxt)
              explain.ifnot(all(iradii >= 0), ctxt)
              explain.ifnot(all(hradii >= 0), ctxt)
              explain.ifnot(all(iradii >= hradii), ctxt)

              par <- list(beta=beta,gamma=gamma,iradii=iradii,hradii=hradii)
              return(par)
            },
            validity=function(par, kind) {
              gamma <- par$gamma
              iradii <- par$iradii
              hradii <- par$hradii
              dh <- diag(hradii)
              dg <- diag(gamma)
              dr <- diag(iradii)
              hhard <- !is.na(dh) & (dh > 0)
              ihard <- !is.na(dr) & (dr > 0) & !is.na(dg) & (dg == 0)
              hard <- hhard | ihard
              operative <- !is.na(gamma) & !is.na(iradii) & (iradii > 0)
              switch(kind,
                     stabilising=all(hard),
                     integrable={
                       all(hard) || all(gamma[operative] <= 1)
                     })
            },
            explainvalid=list(
              integrable=paste(
                "hradii[i,i] > 0 or gamma[i,i] = 0 for all i, or",
                "gamma[i,j] <= 1 for all i and j"),
              stabilising="hradii[i,i] > 0 or gamma[i,i] = 0 for all i"),
            reach=function(par, ...) {
              r <- par$iradii
              h <- par$hradii
              g <- par$gamma
              roperative <- ! (is.na(r) | (g == 1))
              hoperative <- ! is.na(h)
              return(max(0, r[roperative], h[hoperative]))
            }
            ),
#       
# 6. Diggle-Gates-Stibbard interaction
#    (function number 1 from Diggle, Gates, and Stibbard)
       
       'dgs' = 
       list(
            C.id="dgs",
            multitype=FALSE,
            parhandler=function(par, ...) {
              ctxt <- "For the dgs cif"
              par <- check.named.list(par, c("beta","rho"), ctxt)
              with(par, check.finite(beta, ctxt))
              with(par, check.finite(rho, ctxt))
              with(par, explain.ifnot(all(beta >= 0), ctxt))
              with(par, check.1.real(rho, ctxt))
              with(par, explain.ifnot(rho >= 0, ctxt))
              return(par)
            },
            validity=function(par, kind) {
              switch(kind,
                     integrable=TRUE,
                     stabilising=FALSE)
            },
            explainvalid=list(integrable="TRUE", stabilising="FALSE"),
            reach=function(par, ...) {
              return(par[["rho"]])
            }
            ),
#
# 7. Diggle-Gratton interaction 
#    (function number 2 from Diggle, Gates, and Stibbard).

       'diggra' =
       list(
            C.id="diggra",
            multitype=FALSE,
            parhandler=function(par, ...) {
              ctxt <- "For the diggra cif"
              par <- check.named.list(par, c("beta","kappa","delta","rho"),
                                      ctxt)
              with(par, check.finite(beta, ctxt))
              with(par, check.finite(kappa, ctxt))
              with(par, check.finite(delta, ctxt))
              with(par, check.finite(rho, ctxt))
              with(par, explain.ifnot(all(beta >= 0), ctxt))
              with(par, check.1.real(kappa, ctxt))
              with(par, check.1.real(delta, ctxt))
              with(par, check.1.real(rho,   ctxt))
              with(par, explain.ifnot(kappa >= 0, ctxt))              
              with(par, explain.ifnot(delta >= 0, ctxt))              
              with(par, explain.ifnot(rho >= 0, ctxt))              
              with(par, explain.ifnot(delta < rho, ctxt))              
              return(par)
            },
            validity=function(par, kind) {
              switch(kind,
                     integrable=TRUE,
                     stabilising=FALSE)
            },
            explainvalid=list(integrable="TRUE",stabilising="FALSE"),
            reach=function(par, ...) {
              return(par[["rho"]])
            }
            ),
#       
# 8. Geyer saturation model
#       
       'geyer' = 
       list(
            C.id="geyer",
            multitype=FALSE,
            parhandler=function(par, ...) {
              ctxt <- "For the geyer cif"
              par <- check.named.list(par, c("beta","gamma","r","sat"), ctxt)
              with(par, check.1.real(gamma, ctxt))
              with(par, check.1.real(r,     ctxt))
              with(par, check.1.real(sat,   ctxt))
              par <- within(par, sat <- min(sat, .Machine$integer.max-100))
              par <- within(par, if(is.na(gamma)) { r <- 0; gamma <- 1 })
              with(par, check.finite(beta, ctxt))
              with(par, check.finite(gamma, ctxt))
              with(par, check.finite(r, ctxt))
              with(par, check.finite(sat, ctxt))
              with(par, explain.ifnot(all(beta >= 0), ctxt))
              return(par)
            },
            validity=function(par, kind) {
              switch(kind,
                     integrable=TRUE,
                     stabilising=FALSE)
            },
            explainvalid=list(integrable="TRUE", stabilising="FALSE"),
            reach = function(par, ...) {
              r <- par[["r"]]
              g <- par[["gamma"]]
              return(if(g == 1) 0 else 2 * r)
            }
            ),
#       
# 9. The ``lookup'' device.  This permits simulating, at least
# approximately, ANY pairwise interaction function model
# with isotropic pair interaction (i.e. depending only on distance).
# The pair interaction function is provided as a vector of
# distances and corresponding function values which are used
# as a ``lookup table'' by the C code.
#
       'lookup' = 
       list(
            C.id="lookup",
            multitype=FALSE,
            parhandler=function(par, ...) {
              ctxt <- "For the lookup cif"
              par <- check.named.list(par, c("beta","h"), ctxt, "r")
              with(par, check.finite(beta, ctxt))
              with(par, explain.ifnot(all(beta >= 0), ctxt))
              beta   <- par[["beta"]]
              h.init <- par[["h"]]
              r      <- par[["r"]]
              if(is.null(r)) {
		if(!is.stepfun(h.init))
                  stop(paste("For cif=lookup, if component r of",
                             "par is absent then component h must",
                             "be a stepfun object."))
		if(!is.cadlag(h.init))
                  stop(paste("The lookup pairwise interaction step",
			     "function must be right continuous,\n",
			     "i.e. built using the default values of the",
                             sQuote("f"), "and", sQuote("right"),
                             "arguments for stepfun."))
		r     <- knots(h.init)
		h0    <- get("yleft",envir=environment(h.init))
		h     <- h.init(r)
		nlook <- length(r)
		if(!identical(all.equal(h[nlook],1),TRUE))
                  stop(paste("The lookup interaction step function",
                             "must be equal to 1 for", dQuote("large"),
                             "distances."))
		if(r[1] <= 0)
                  stop(paste("The first jump point (knot) of the lookup",
                             "interaction step function must be",
                             "strictly positive."))
		h <- c(h0,h)
              } else {
		h     <- h.init
		nlook <- length(r)
		if(length(h) != nlook)
                  stop("Mismatch of lengths of h and r lookup vectors.")
		if(any(is.na(r)))
                  stop("Missing values not allowed in r lookup vector.")
		if(is.unsorted(r))
                  stop("The r lookup vector must be in increasing order.")
		if(r[1] <= 0)
                  stop(paste("The first entry of the lookup vector r",
                             "should be strictly positive."))
		h <- c(h,1)
              }
              if(any(h < 0))
		stop(paste("Negative values in the lookup",
                           "pairwise interaction function."))
              if(h[1] > 0 & any(h > 1))
		stop(paste("Lookup pairwise interaction function does",
                           "not define a valid point process."))
              rmax   <- r[nlook]
              r <- c(0,r)
              nlook <- nlook+1
              deltar <- mean(diff(r))
              if(identical(all.equal(diff(r),rep.int(deltar,nlook-1)),TRUE)) {
		par <- list(beta=beta,nlook=nlook,
                            equisp=1,
                            deltar=deltar,rmax=rmax, h=h)
              } else {
		par <- list(beta=beta,nlook=nlook,
                            equisp=0,
                            deltar=deltar,rmax=rmax, h=h,
                            r=r)
              }
              return(par) 
            },
            validity=function(par, kind) {
              h <- par$h
              if(is.stepfun(h))
                h <- eval(expression(c(yleft,y)),envir=environment(h))
              switch(kind,
                     integrable={
                       (h[1] == 0) || all(h <= 1)
                     },
                     stabilising={ h[1] == 0 })
            },
            explainvalid=list(
              integrable="h[1] == 0 or h[i] <= 1 for all i",
              stabilising="h[1] == 0"),
            reach = function(par, ...) {
              r <- par[["r"]]
              h <- par[["h"]]
              if(is.null(r)) 
                r <- knots(h)
              return(max(r))
            }
            ),
#       
# 10. Area interaction
#       
       'areaint'=
       list(
            C.id="areaint",
            multitype=FALSE,
            parhandler=function(par, ...) {
              ctxt <- "For the areaint cif"
              par <- check.named.list(par, c("beta","eta","r"), ctxt)
              par <- within(par, if(is.na(r)) { r <- 0 })
              with(par, check.finite(beta, ctxt))
              with(par, explain.ifnot(all(beta >= 0), ctxt))
              with(par, check.1.real(eta, ctxt))
              with(par, check.1.real(r,   ctxt))
              with(par, check.finite(eta, ctxt))
              with(par, check.finite(r,   ctxt))
              with(par, explain.ifnot(eta >= 0, ctxt))
              with(par, explain.ifnot(r >= 0,   ctxt))
              return(par)
            },
            validity=function(par, kind) {
              switch(kind,
                     integrable=TRUE,
                     stabilising=FALSE)
            },
            explainvalid=list(integrable="TRUE", stabilising="FALSE"),
            reach = function(par, ...) {
              r <- par[["r"]]
              eta <- par[["eta"]]
              return(if(eta == 1) 0 else (2 * r))
            }
            ),
#
# 11. The ``badgey'' (Baddeley-Geyer) model.
#
       'badgey' =
       list(
            C.id="badgey",
            multitype=FALSE,
            parhandler=function(par, ...) {
              ctxt <- "For the badgey cif"
              par <- check.named.list(par, c("beta","gamma","r","sat"), ctxt)
              par <- within(par, sat <- pmin(sat, .Machine$integer.max-100))
              par <- within(par, gamma[is.na(gamma) | is.na(r)] <- 1)
              par <- within(par, r[is.na(r)] <- 0)
              with(par, check.finite(beta, ctxt))
              with(par, check.finite(gamma, ctxt))
              with(par, check.finite(r, ctxt))
              with(par, check.finite(sat, ctxt))
              with(par, explain.ifnot(all(beta >= 0), ctxt))
              with(par, explain.ifnot(all(gamma >= 0), ctxt))
              with(par, explain.ifnot(all(r >= 0), ctxt))
              with(par, explain.ifnot(all(sat >= 0), ctxt))
              with(par, explain.ifnot(length(gamma) == length(r), ctxt)) 
              gamma <- par[["gamma"]]
              r     <- par[["r"]]
              sat   <- par[["sat"]]
              if(length(sat)==1) sat <- rep.int(sat,length(gamma))
              else explain.ifnot(length(sat) == length(gamma), ctxt)
              mmm <- cbind(gamma,r,sat)
              mmm <- mmm[fave.order(r),]
              ndisc <- length(r)
              par <- list(beta=par$beta,ndisc=ndisc,parms=as.vector(t(mmm)))
              return(par)
            },
            validity=function(par, kind) {
              switch(kind,
                     integrable=TRUE,
                     stabilising=FALSE)
            },
            explainvalid=list(integrable="TRUE", stabilising="FALSE"),
            reach = function(par, ...) {
              r <- par[["r"]]
              gamma <- par[["gamma"]]
              operative <- (gamma != 1)
              return(if(!any(operative)) 0 else (2 * max(r[operative])))
            }
            ),
#
# 12. The hard core process
       'hardcore' =
       list(
            C.id="hardcore",
            multitype=FALSE,
            parhandler=function(par, ...) {
              ctxt <- "For the hardcore cif"
              par <- check.named.list(par, c("beta", "hc"), ctxt)
              par <- within(par, if(is.na(hc)) { hc <- 0 })
              with(par, check.finite(beta, ctxt))
              with(par, check.finite(hc, ctxt))
              with(par, explain.ifnot(all(beta >= 0), ctxt))
              with(par, check.1.real(hc, ctxt))
              return(par)
            },
            validity=function(par, kind) {
              hc <- par$hc
              switch(kind,
                     integrable=TRUE,
                     stabilising=(hc > 0))
            },
            explainvalid=list(integrable="TRUE", stabilising="hc > 0"),
            reach = function(par, ...) {
              hc <- par[["hc"]]
              return(hc)
            }
            ),
#
# Lucky 13. Fiksel process
       'fiksel' =
       list(
            C.id="fiksel",
            multitype=FALSE,
            parhandler=function(par, ...) {
              ctxt <- "For the Fiksel cif"
              par <- check.named.list(par,
                                      c("beta", "r", "hc", "kappa", "a"),
                                      ctxt)
              with(par, check.finite(beta, ctxt))
              with(par, check.finite(r, ctxt))
              with(par, check.finite(hc, ctxt))
              with(par, check.finite(kappa, ctxt))
              with(par, check.finite(a, ctxt))
              with(par, check.1.real(r, ctxt))
              with(par, check.1.real(hc, ctxt))
              with(par, check.1.real(kappa, ctxt))
              with(par, check.1.real(a, ctxt))
              with(par, explain.ifnot(all(beta >= 0), ctxt))
              with(par, explain.ifnot(hc >= 0, ctxt))
              with(par, explain.ifnot(r > hc, ctxt))
              return(par)
            },
            validity=function(par, kind) {
              hc <- par$hc
              a  <- par$a
              switch(kind,
                     integrable=(hc > 0 || a <= 0),
                     stabilising=(hc > 0))
            },
            explainvalid=list(
              integrable="hc > 0 or a <= 0",
              stabilising="hc > 0"),
            reach = function(par, ...) {
              r <- par[["r"]]
              hc <- par[["hc"]]              
              a <- par[["a"]]
              return(if(a != 0) r else hc)
            }
            ),
#
# 14. Lennard-Jones
       'lennard' =
       list(
            C.id="lennard",
            multitype=FALSE,
            parhandler=function(par, ...) {
              ctxt <- "For the Lennard-Jones cif"
              par <- check.named.list(par,
                                      c("beta", "sigma", "epsilon"),
                                      ctxt)
              with(par, check.finite(beta, ctxt))
              with(par, check.finite(sigma, ctxt))
              with(par, check.finite(epsilon, ctxt))
              with(par, explain.ifnot(all(beta >= 0), ctxt))
              with(par, check.1.real(sigma, ctxt))
              with(par, check.1.real(epsilon, ctxt))
              with(par, explain.ifnot(sigma > 0, ctxt))
              with(par, explain.ifnot(epsilon > 0, ctxt))
              return(par)
            },
            validity=function(par, kind) {
              switch(kind,
                     integrable=(par$sigma > 0),
                     stabilising=FALSE)
            },
            explainvalid=list(
              integrable="sigma > 0",
              stabilising="FALSE"),
            reach = function(par, ...) {
              sigma <- par[["sigma"]]
              return(2.5 * sigma)
            }
            ),
#       
# 15. Multitype hardcore.
#       
       'multihard' = 
       list(
            C.id="multihard",
            multitype=TRUE,
            parhandler=function(par, types) {
              ctxt="For the multihard cif"
              par <- check.named.list(par,
                                      c("beta","hradii"),
                                      ctxt)

              beta <- par$beta
              hradii <- par$hradii
              ntypes <- length(types)

              check.nvector(beta, ntypes, TRUE, "types")
              check.finite(beta, ctxt)
              
              MultiPair.checkmatrix(hradii, ntypes, "par$hradii")
              hradii[is.na(hradii)] <- 0
              check.finite(hradii, ctxt)

              explain.ifnot(all(beta >= 0), ctxt)
              explain.ifnot(all(hradii >= 0), ctxt)

              par <- list(beta=beta,hradii=hradii)
              return(par)
            },
            validity=function(par, kind) {
              switch(kind,
                     integrable=return(TRUE),
                     stabilising={
                       hself <- diag(par$hradii)
                       repel <- !is.na(hself) & (hself > 0)
                       return(all(repel))
                     })
            },
            explainvalid=list(
              integrable="TRUE",
              stabilising="hradii[i,i] > 0 for all i"),
            reach=function(par, ...) {
              return(max(0, par$hradii, na.rm=TRUE))
            }
            ),
#       
# 16. Triplets.
#       
       'triplets'=
       list(
            C.id="triplets",
            multitype=FALSE,
            parhandler=function(par, ...) {
              ctxt <- "For the triplets cif"
              par <- check.named.list(par, c("beta","gamma","r"), ctxt)
              # treat r=NA as absence of interaction
              par <- within(par, if(is.na(r)) { r <- 0; gamma <- 1 })
              with(par, check.finite(beta, ctxt))
              with(par, check.finite(gamma, ctxt))
              with(par, check.finite(r, ctxt))
              with(par, check.1.real(gamma, ctxt))
              with(par, check.1.real(r,     ctxt))
              with(par, explain.ifnot(all(beta >= 0), ctxt))
              with(par, explain.ifnot(gamma >= 0, ctxt))
              with(par, explain.ifnot(r >= 0, ctxt))
              return(par)
            },
            validity=function(par, kind) {
              gamma <- par$gamma
              switch(kind,
                     integrable=(gamma <= 1),
                     stabilising=(gamma == 0)
                     )
            },
            explainvalid=list(
              integrable="gamma <= 1",
              stabilising="gamma == 0"),
            reach = function(par, ...) {
              r <- par[["r"]]
              g <- par[["gamma"]]
              return(if(g == 1) 0 else r)
            }
            )
       # end of list '.Spatstat.RmhTable'
       )


