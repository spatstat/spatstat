#
#  rmhmodel.ppm.R
#
#   convert ppm object into format palatable to rmh.default
#
#  $Revision: 2.58 $   $Date: 2014/09/16 09:13:55 $
#
#   .Spatstat.rmhinfo
#   rmhmodel.ppm()
#

.Spatstat.Rmhinfo <-
list(
     "Multitype Hardcore process" =
     function(coeffs, inte) {
       # hard core radii r[i,j]
       hradii <- inte$par$hradii
       return(list(cif='multihard',
                   par=list(hradii=hradii),
                   ntypes=ncol(hradii)))
     },
     "Lennard-Jones process" =
     function(coeffs, inte) {
       sigma   <- inte$par$sigma
       epsilon <- inte$par$epsilon
       return(list(cif='lennard',
                   par=list(sigma=sigma, epsilon=epsilon),
                   ntypes=1))
     },
     "Fiksel process" =
     function(coeffs, inte) {
       hc <- inte$par$hc
       r  <- inte$par$r
       kappa <- inte$par$kappa
       a <- inte$interpret(coeffs,inte)$param$a
       return(list(cif='fiksel',
                   par=list(r=r,hc=hc,kappa=kappa,a=a),
                   ntypes=1))
     },
     "Diggle-Gates-Stibbard process" =
     function(coeffs, inte) {
       rho   <- inte$par$rho
       return(list(cif='dgs',
                   par=list(rho=rho),
                   ntypes=1))
     },
     "Diggle-Gratton process" =
     function(coeffs, inte) {
       kappa <- inte$interpret(coeffs,inte)$param$kappa
       delta <- inte$par$delta
       rho   <- inte$par$rho
       return(list(cif='diggra',
                   par=list(kappa=kappa,delta=delta,rho=rho),
                   ntypes=1))
     },
     "Hard core process" =
     function(coeffs, inte) {
       hc <- inte$par$hc
       return(list(cif='hardcore',
                   par=list(hc=hc),
                   ntypes=1))
     },
     "Geyer saturation process" =
     function(coeffs, inte) {
       gamma <- inte$interpret(coeffs,inte)$param$gamma
       r <- inte$par$r
       sat <- inte$par$sat
       return(list(cif='geyer',
                   par=list(gamma=gamma,r=r,sat=sat),
                   ntypes=1))
     },
     "Soft core process" =
     function(coeffs, inte) {
       kappa <- inte$par$kappa
       sigma <- inte$interpret(coeffs,inte)$param$sigma
       return(list(cif="sftcr",
                   par=list(sigma=sigma,kappa=kappa),
                   ntypes=1))
     },
     "Strauss process" =
     function(coeffs, inte) {
       gamma <- inte$interpret(coeffs,inte)$param$gamma
       r <- inte$par$r
       return(list(cif = "strauss",
                   par = list(gamma = gamma, r = r),
                   ntypes=1))
     },
     "Strauss - hard core process" =
     function(coeffs, inte) {
       gamma <- inte$interpret(coeffs,inte)$param$gamma
       r <- inte$par$r
       hc <- inte$par$hc
       return(list(cif='straush',
                   par=list(gamma=gamma,r=r,hc=hc),
                   ntypes=1))
     },
     "Triplets process" =
     function(coeffs, inte) {
       gamma <- inte$interpret(coeffs,inte)$param$gamma
       r <- inte$par$r
       return(list(cif = "triplets",
                   par = list(gamma = gamma, r = r),
                   ntypes=1))
     },
     "Multitype Strauss process" =
     function(coeffs, inte) {
       # interaction radii r[i,j]
       radii <- inte$par$radii
       # interaction parameters gamma[i,j]
       gamma <- (inte$interpret)(coeffs, inte)$param$gammas
       return(list(cif='straussm',
                   par=list(gamma=gamma,radii=radii),
                   ntypes=ncol(radii)))
     },
     "Multitype Strauss Hardcore process" =
     function(coeffs, inte) {
       # interaction radii r[i,j]
       iradii <- inte$par$iradii
       # hard core radii r[i,j]
       hradii <- inte$par$hradii
       # interaction parameters gamma[i,j]
       gamma <- (inte$interpret)(coeffs, inte)$param$gammas
       return(list(cif='straushm',
                   par=list(gamma=gamma,iradii=iradii,hradii=hradii),
                   ntypes=ncol(iradii)))
     },
     "Piecewise constant pairwise interaction process" =
     function(coeffs, inte) {
       r <- inte$par$r
       gamma <- (inte$interpret)(coeffs, inte)$param$gammas
       h <- stepfun(r, c(gamma, 1))
       return(list(cif='lookup', par=list(h=h),
                   ntypes=1))
     },
     "Area-interaction process" =
     function(coeffs, inte) {
       r <- inte$par$r
       eta <- (inte$interpret)(coeffs, inte)$param$eta
       return(list(cif='areaint', par=list(eta=eta,r=r), ntypes=1))
     },
     "hybrid Geyer process" =
     function(coeffs, inte) {
       r <- inte$par$r
       sat <- inte$par$sat
       gamma <- (inte$interpret)(coeffs,inte)$param$gammas
       return(list(cif='badgey',par=list(gamma=gamma,r=r,sat=sat), ntypes=1))
     },
     "Hybrid interaction"=
     function(coeffs, inte){
       # for hybrids, $par is a list of the component interactions
       interlist <- inte$par
       # check for Poisson components
       ispois <- unlist(lapply(interlist, is.poisson))
       if(all(ispois)) {
         # reduces to Poisson
         Z <- list(cif='poisson', par=list())
         return(Z)
       } else if(any(ispois)) {
         # remove Poisson components
         interlist <- interlist[!ispois]
       }
       # 
       N <- length(interlist)
       cifs <- character(N)
       pars <- vector(mode="list", length=N)
       ntyp <- integer(N)
       for(i in 1:N) {
         interI <- interlist[[i]]
         # forbid hybrids-of-hybrids - these should not occur anyway
         if(interI$name == "Hybrid interaction")
           stop("Simulation of a hybrid-of-hybrid interaction is not implemented")
         # get RMH mapping for I-th component
         siminfoI <- .Spatstat.Rmhinfo[[interI$name]]
         if(is.null(siminfoI))
           stop(paste("Simulation of a fitted", sQuote(interI$name),
                      "has not yet been implemented"),
                call.=FALSE)
         # nameI is the tag that identifies I-th component in hybrid
         nameI  <- names(interlist)[[i]]
         nameI. <- paste(nameI, ".", sep="")
         # find coefficients with prefix that exactly matches nameI.
         Cname  <- names(coeffs)
         prefixlength <- nchar(nameI.)
         Cprefix <- substr(Cname, 1, prefixlength)
         relevant <- (Cprefix == nameI.)
         # extract coefficients
         #   (there may be none, if this interaction is an 'offset')
         coeffsI <- coeffs[relevant]
         # remove the prefix so the coefficients are recognisable to 'siminfoI'
         if(any(relevant)) 
           names(coeffsI) <-
             substr(Cname[relevant], prefixlength+1, max(nchar(Cname)))
         # compute RMH info
         ZI <- siminfoI(coeffsI, interI)
         cifs[i] <- ZI$cif
         pars[[i]] <- ZI$par
         ntyp[i] <- ZI$ntypes
       }
       nt <- unique(ntyp[ntyp != 1])
       if(length(nt) > 1)
         stop(paste("Hybrid components have different numbers of types:",
                    commasep(nt)))
       if(N == 1) {
         # single cif: revert to original format: par is a list of parameters
         Z <- list(cif=cifs[1], par=pars[[1]], ntypes=ntyp)
       } else {
         # hybrid cif: par is a list of lists of parameters
         Z <- list(cif=cifs,    par=pars,      ntypes=ntyp)
       }
       return(Z)
     }
)


# OTHER MODELS not yet implemented:
#
#
#      interaction object           rmh.default 
#      ------------------           -----------
#
#           OrdThresh                <none>
#


rmhmodel.ppm <- function(model, win, ...,
                         verbose=TRUE, project=TRUE,
                         control=rmhcontrol(),
                         new.coef=NULL) {
  ## converts ppm object `model' into format palatable to rmh.default
  
  verifyclass(model, "ppm")

  if(!is.null(new.coef)) {
    ## hack the coefficients
    co <- coef(model)
    check.nvector(new.coef, length(co), things="coefficients")
    model$coef.orig <- co
    model$coef <- new.coef
  }

  ## Ensure the fitted model is valid
  ## (i.e. exists mathematically as a point process)
  if(!valid.ppm(model)) {
    if(project) {
      if(verbose)
        cat("Model is invalid - projecting it\n")
      model <- project.ppm(model, fatal=TRUE)
    } else stop("The fitted model is not a valid point process")
  }
    
  if(verbose)
    cat("Extracting model information...")
    
  ## Extract essential information
  Y <- summary(model, quick="no variances")

  if(Y$marked && !Y$multitype)
    stop("Not implemented for marked point processes other than multitype")

  if(Y$uses.covars && is.data.frame(model$covariates))
    stop(paste("This model cannot be simulated, because the",
               "covariate values were given as a data frame."))
    
  ## enforce defaults for `control'

  control <- rmhcontrol(control)

  ## adjust to peculiarities of model
    
  control <- rmhResolveControl(control, model)
    
  ########  Interpoint interaction
  if(Y$poisson) {
    Z <- list(cif="poisson",
              par=list())  # par is filled in later
  } else {
    ## First check version number of ppm object
    if(Y$antiquated) 
      stop(paste("This model was fitted by a very old version",
                 "of the package: spatstat", Y$version,
                 "; simulation is not possible.",
                 "Re-fit the model using your original code"))
    else if(Y$old)
      warning(paste("This model was fitted by an old version",
                    "of the package: spatstat", Y$version,
                    ". Re-fit the model using update.ppm",
                    "or your original code"))
    ## Extract the interpoint interaction object
    inte <- Y$entries$interaction
    ## Determine whether the model can be simulated using rmh
    siminfo <- .Spatstat.Rmhinfo[[inte$name]]
    if(is.null(siminfo))
      stop(paste("Simulation of a fitted", sQuote(inte$name),
                 "has not yet been implemented"))
      
    ## Get fitted model's canonical coefficients
    coeffs <- Y$entries$coef
    if(newstyle.coeff.handling(inte)) {
      ## extract only the interaction coefficients
      Vnames <- Y$entries$Vnames
      IsOffset <- Y$entries$IsOffset
      coeffs <- coeffs[Vnames[!IsOffset]]
    }
    ## Translate the model to the format required by rmh.default
    Z <- siminfo(coeffs, inte)
    if(is.null(Z))
      stop("The model cannot be simulated")
    else if(is.null(Z$cif))
      stop(paste("Internal error: no cif returned from .Spatstat.Rmhinfo"))
  }

  ## Don't forget the types
  if(Y$multitype && is.null(Z$types))
    Z$types <- levels(Y$entries$marks)
       
  ######## Window for result 
    
  if(missing(win))
    win <- Y$entries$data$window

  Z$w <- win

  ######## Expanded window for simulation?

  covims <- if(Y$uses.covars) model$covariates[Y$covars.used] else NULL
    
  wsim <- rmhResolveExpansion(win, control, covims, "covariate")$wsim
      
  ###### Trend or Intensity ############

  if(verbose)
    cat("Evaluating trend...")
    
  if(Y$stationary) {
    ## first order terms (beta or beta[i]) are carried in Z$par
    beta <- as.numeric(Y$trend$value)
    Z$trend <- NULL
  } else {
    ## trend terms present
    ## all first order effects are subsumed in Z$trend
    beta <- if(!Y$marked) 1 else rep.int(1, length(Z$types))
    ## predict on window possibly larger than original data window
    Z$trend <- 
      if(wsim$type == "mask")
        predict(model, window=wsim, type="trend", locations=wsim)
      else 
        predict(model, window=wsim, type="trend")
  }
    
  Ncif <- length(Z$cif)
  if(Ncif == 1) {
    ## single interaction
    Z$par[["beta"]] <- beta
  } else {
    ## hybrid interaction
    if(all(Z$ntypes == 1)) {
      ## unmarked model: scalar 'beta' is absorbed in first cif
      absorb <- 1
    } else {
      ## multitype model: vector 'beta' is absorbed in a multitype cif
      absorb <- min(which(Z$ntypes > 1))
    }
    Z$par[[absorb]]$beta <- beta
    ## other cifs have par$beta = 1 
    for(i in (1:Ncif)[-absorb])
      Z$par[[i]]$beta <- rep.int(1, Z$ntypes[i])
  }
    
  if(verbose)
    cat("done.\n")
  Z <- rmhmodel(Z, ...)
  return(Z)
}

rmhResolveExpansion <- function(win, control, imagelist, itype="covariate") {
  # Determine expansion window for simulation
  ex <- control$expand
  
# The following is redundant because it is implied by !will.expand(ex)  
#  if(ex$force.noexp) {
#    # Expansion prohibited
#    return(list(wsim=win, expanded=FALSE))
#  }
  
  # Is expansion contemplated?
  if(!will.expand(ex))
    return(list(wsim=win, expanded=FALSE))

  # Proposed expansion window
  wexp <- expand.owin(win, ex)

  # Check feasibility
  isim <- unlist(lapply(imagelist, is.im))
  imagelist <- imagelist[isim]

  if(length(imagelist) == 0) {
    # Unlimited expansion is feasible
    return(list(wsim=wexp, expanded=TRUE))
  }

  # Expansion is limited to domain of image data
  # Determine maximum possible expansion window
  wins <- lapply(imagelist, as.owin)
  cwin <- do.call("intersect.owin", unname(wins))
  
  if(!is.subset.owin(wexp, cwin)) {
    # Cannot expand to proposed window
    if(ex$force.exp)
      stop(paste("Cannot expand the simulation window,",
                 "because the", itype, "images do not cover",
                 "the expanded window"), call.=FALSE)
      # Take largest possible window
    wexp <- intersect.owin(wexp, cwin)
  }
  return(list(wsim=wexp, expanded=TRUE))
}

