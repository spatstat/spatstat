#    mpl.R
#
#	$Revision: 5.233 $	$Date: 2020/03/04 04:25:13 $
#
#    mpl.engine()
#          Fit a point process model to a two-dimensional point pattern
#          by maximum pseudolikelihood
#
#    mpl.prepare()
#          set up data for glm procedure
#
# -------------------------------------------------------------------
#

# "mpl" <- function(Q,
#         trend = ~1,
#	 interaction = NULL,
#         data = NULL,
#	 correction="border",
#	 rbord = 0,
#         use.gam=FALSE) {
#   .Deprecated("ppm", package="spatstat")
#   ppm(Q=Q, trend=trend, interaction=interaction,
#       covariates=data, correction=correction, rbord=rbord,
#       use.gam=use.gam, method="mpl")
# }

mpl.engine <- 
  function(Q,
           trend = ~1,
           interaction = NULL,
           ...,
           covariates = NULL,
           subsetexpr = NULL,
           clipwin = NULL,
           covfunargs = list(),
           correction="border",
           rbord = 0,
           use.gam=FALSE,
           gcontrol=list(),
           GLM=NULL,
           GLMfamily=NULL,
           GLMcontrol=NULL,
           famille=NULL,
           forcefit=FALSE,
           nd = NULL,
           eps = eps,
           allcovar=FALSE,
           callstring="",
           precomputed=NULL,
           savecomputed=FALSE,
           preponly=FALSE,
           rename.intercept=TRUE,
           justQ = FALSE,
           weightfactor = NULL)
  {
    GLMname <- if(!missing(GLM)) short.deparse(substitute(GLM)) else NULL
    ## Extract precomputed data if available
    if(!is.null(precomputed$Q)) {
      Q <- precomputed$Q
      X <- precomputed$X
      P <- precomputed$U
    } else {
      ## Determine quadrature scheme from argument Q
      if(verifyclass(Q, "quad", fatal=FALSE)) {
        ## user-supplied quadrature scheme - validate it
        validate.quad(Q, fatal=TRUE, repair=FALSE, announce=TRUE)
        ## Extract data points
        X <- Q$data
      } else if(verifyclass(Q, "ppp", fatal = FALSE)) {
        ## point pattern - create default quadrature scheme
        X <- Q
        Q <- quadscheme(X, nd=nd, eps=eps, check=FALSE)
      } else 
      stop("First argument Q should be a point pattern or a quadrature scheme")
      ## Data and dummy points together
      P <- union.quad(Q)
    }
    ## clip to subset?
    if(!is.null(clipwin)) {
      if(is.data.frame(covariates)) 
        covariates <- covariates[inside.owin(P, w=clipwin), , drop=FALSE]
      Q <- Q[clipwin]
      X <- X[clipwin]
      P <- P[clipwin]
    }
    ## secret exit  
    if(justQ) return(Q)
    ##  
    computed <- if(savecomputed) list(X=X, Q=Q, U=P) else NULL
    ##
    ## Validate main arguments
    if(!is.null(trend) && !inherits(trend, "formula"))
      stop(paste("Argument", sQuote("trend"), "must be a formula"))
    if(!is.null(interaction) && !inherits(interaction, "interact"))
      stop(paste("Argument", sQuote("interaction"), "has incorrect format"))
    ##
    check.1.real(rbord, "In ppm")
    explain.ifnot(rbord >= 0, "In ppm")
    ## rbord applies only to border correction
    if(correction != "border") rbord <- 0 
    ##
    covfunargs <- as.list(covfunargs)
    ##
    ## Interpret the call
    if(is.null(trend)) {
      trend <- ~1
      environment(trend) <- parent.frame()
    }
    want.trend <- !identical.formulae(trend, ~1)
    want.inter <- !is.null(interaction) && !is.null(interaction$family)

    ## Stamp with spatstat version number
    spv <- package_version(versionstring.spatstat())
    the.version <- list(major=spv$major,
                        minor=spv$minor,
                        release=spv$patchlevel,
                        date="$Date: 2020/03/04 04:25:13 $")

    if(want.inter) {
      ## ensure we're using the latest version of the interaction object
      if(outdated.interact(interaction)) 
        interaction <- update(interaction)
    }

    ##  
  
    if(!want.trend && !want.inter &&
       !forcefit && !allcovar && is.null(subsetexpr)) {
      ## the model is the uniform Poisson process
      ## The MPLE (= MLE) can be evaluated directly
      npts <- npoints(X)
      W    <- as.owin(X)
      if(correction == "border" && rbord > 0) {
        npts <- sum(bdist.points(X) >= rbord)
        areaW <- eroded.areas(W, rbord)
      } else {
        npts <- npoints(X)
        areaW <- area(W)
      }
      volume <- areaW * markspace.integral(X)
      lambda <- npts/volume
      ## fitted canonical coefficient
      co <- log(lambda)
      ## asymptotic variance of canonical coefficient
      varcov <- matrix(1/npts, 1, 1)
      fisher <- matrix(npts,   1, 1)
      se <- sqrt(1/npts)
      ## give names
      tag <- if(rename.intercept) "log(lambda)" else "(Intercept)"
      names(co) <- tag
      dimnames(varcov) <- dimnames(fisher) <- list(tag, tag)
      ## maximised log likelihood
      maxlogpl <- if(npts == 0) 0 else npts * (log(lambda) - 1)
      ##
      rslt <- list(
                   method      = "mpl",
                   fitter      = "exact",
                   projected   = FALSE,
                   coef        = co,
                   trend       = trend,
                   interaction = NULL,
                   fitin       = fii(),
                   Q           = Q,
                   maxlogpl    = maxlogpl,
                   satlogpl    = NULL,
                   internal    = list(computed=computed, se=se),
                   covariates  = mpl.usable(covariates),
                                        ## covariates are still retained!
                   covfunargs  = covfunargs,
                   subsetexpr  = NULL,
                   correction  = correction,
                   rbord       = rbord,
                   terms       = terms(trend),
                   fisher      = fisher,
                   varcov      = varcov,
                   version     = the.version,
                   problems    = list())
      class(rslt) <- "ppm"
      return(rslt)
    }
    #################  P r e p a r e    D a t a   ######################

    prep <- mpl.prepare(Q, X, P, trend, interaction,
                        covariates, 
                        want.trend, want.inter, correction, rbord,
                        "quadrature points", callstring,
                        subsetexpr=subsetexpr,
                        allcovar=allcovar,
                        precomputed=precomputed, savecomputed=savecomputed,
                        covfunargs=covfunargs,
                        weightfactor=weightfactor,
                        ...)
    ## back door
    if(preponly) {
      ## exit now, returning prepared data frame and internal information
      prep$info <- list(want.trend=want.trend,
                        want.inter=want.inter,
                        correction=correction,
                        rbord=rbord,
                        interaction=interaction)
      return(prep)
    }
  
  
    fmla <- prep$fmla
    glmdata <- prep$glmdata
    problems <- prep$problems
    likelihood.is.zero <- prep$likelihood.is.zero
    is.identifiable <- prep$is.identifiable
    computed <- resolve.defaults(prep$computed, computed)
    IsOffset <- prep$IsOffset

    ## update covariates (if they were resolved from the environment)  
    if(!is.null(prep$covariates))
      covariates <- prep$covariates
  
    ################# F i t    i t   ####################################

    if(!is.identifiable) 
      stop(paste("in", callstring, ":", problems$unidentifiable$print),
           call.=FALSE)
  
    ## to avoid problem with package checker  
    .mpl.W <- glmdata$.mpl.W
    .mpl.SUBSET <- glmdata$.mpl.SUBSET

    ## determine algorithm control parameters
    if(is.null(gcontrol)) gcontrol <- list() else stopifnot(is.list(gcontrol))
    gcontrol <- if(!is.null(GLMcontrol)) do.call(GLMcontrol, gcontrol) else
                if(want.trend && use.gam) do.call(mgcv::gam.control, gcontrol) else
                do.call(stats::glm.control, gcontrol)
  
    ## Fit the generalized linear/additive model.

    if(is.null(GLM) && is.null(famille)) {
      ## the sanctioned technique, using `quasi' family
      if(want.trend && use.gam) {
        FIT  <- gam(fmla, family=quasi(link="log", variance="mu"),
                    weights=.mpl.W,
                    data=glmdata, subset=.mpl.SUBSET,
                    control=gcontrol)
        fittername <- "gam"
      } else {
        FIT  <- glm(fmla, family=quasi(link="log", variance="mu"),
                    weights=.mpl.W,
                    data=glmdata, subset=.mpl.SUBSET,
                    control=gcontrol, model=FALSE)
        fittername <- "glm"
      }
    } else if(!is.null(GLM)) {
      ## alternative GLM fitting function or penalised GLM etc
      fam <- GLMfamily %orifnull% quasi(link="log", variance="mu")
      FIT <- GLM(fmla, family=fam,
                 weights=.mpl.W,
                 data=glmdata, subset=.mpl.SUBSET,
                 control=gcontrol)
      fittername <- GLMname
    } else {
      ## experimentation only!
      if(is.function(famille))
        famille <- famille()
      stopifnot(inherits(famille, "family"))
      if(want.trend && use.gam) {
        FIT  <- gam(fmla, family=famille, weights=.mpl.W,
                    data=glmdata, subset=.mpl.SUBSET,
                    control=gcontrol)
        fittername <- "experimental"
      } else {
        FIT  <- glm(fmla, family=famille, weights=.mpl.W,
                    data=glmdata, subset=.mpl.SUBSET,
                    control=gcontrol, model=FALSE)
        fittername <- "experimental"
      }
    }
    environment(FIT$terms) <- sys.frame(sys.nframe())

  
    ################  I n t e r p r e t    f i t   #######################

    ## Fitted coefficients
    co <- FIT$coef

    ## glm covariates
    W <- glmdata$.mpl.W
    SUBSET <- glmdata$.mpl.SUBSET        
    Z <- is.data(Q)
    Vnames <- prep$Vnames
    vnameprefix <- prep$vnameprefix

    ## saturated log pseudolikelihood
    satlogpl <- - (sum(log(W[Z & SUBSET])) + sum(Z & SUBSET))
    ## attained value of max log pseudolikelihood
    maxlogpl <- if(likelihood.is.zero) -Inf else (satlogpl - deviance(FIT)/2)

    ## fitted interaction object
    fitin <- if(want.inter) {
               fii(interaction, co, Vnames, IsOffset, vnameprefix)
             } else fii()
    unitname(fitin) <- unitname(X)
    ######################################################################
    ## Clean up & return 

    rslt <-
      list(
           method       = "mpl",
           fitter       = fittername,
           projected    = FALSE,
           coef         = co,
           trend        = trend,
           interaction  = if(want.inter) interaction else NULL,
           fitin        = fitin,
           Q            = Q,
           maxlogpl     = maxlogpl,
           satlogpl     = satlogpl,
           internal     = list(glmfit=FIT, glmdata=glmdata, Vnames=Vnames,
                               IsOffset=IsOffset, fmla=fmla, computed=computed,
                               vnamebase=prep$vnamebase,
                               vnameprefix=prep$vnameprefix),
           covariates   = mpl.usable(covariates),
           covfunargs   = covfunargs,
           subsetexpr   = subsetexpr,
           correction   = correction,
           rbord        = rbord,
           terms        = terms(trend),
           version      = the.version,
           problems     = problems)
    class(rslt) <- "ppm"
    return(rslt)
  }  



##########################################################################
### /////////////////////////////////////////////////////////////////////
##########################################################################

mpl.prepare <- local({

  mpl.prepare <- function(Q, X, P, trend, interaction, covariates, 
                          want.trend, want.inter, correction, rbord,
                          Pname="quadrature points", callstring="",
                          ...,
                          subsetexpr=NULL,
                          covfunargs=list(),
                          allcovar=FALSE,
                          precomputed=NULL, savecomputed=FALSE,
                          vnamebase=c("Interaction", "Interact."),
                          vnameprefix=NULL,
                          warn.illegal=TRUE,
                          warn.unidentifiable=TRUE,
                          weightfactor=NULL,
                          skip.border=FALSE,
                          clip.interaction=TRUE,
                          splitInf=FALSE) {
    ## Q: quadrature scheme
    ## X = data.quad(Q)
    ## P = union.quad(Q)
  
    if(missing(want.trend))
      want.trend <- !is.null(trend) && !identical.formulae(trend, ~1)
    if(missing(want.inter))
      want.inter <- !is.null(interaction) && !is.null(interaction$family)

    want.subset <- !is.null(subsetexpr)
  
    computed <- list()
    problems <- list()
  
    names.precomputed <- names(precomputed)

    likelihood.is.zero <- FALSE
    is.identifiable <- TRUE
  
    if(!missing(vnamebase)) {
      if(length(vnamebase) == 1)
        vnamebase <- rep.int(vnamebase, 2)
      if(!is.character(vnamebase) || length(vnamebase) != 2)
        stop("Internal error: illegal format of vnamebase")
    }
    if(!is.null(vnameprefix)) {
      if(!is.character(vnameprefix) || length(vnameprefix) != 1)
        stop("Internal error: illegal format of vnameprefix")
    }
      
    ################ C o m p u t e     d a t a  ####################

    ## Extract covariate values
    updatecovariates <- FALSE
    covariates.df <- NULL
    if(allcovar || want.trend || want.subset) {
      if("covariates.df" %in% names.precomputed) {
        covariates.df <- precomputed$covariates.df
      } else {
        if(!is.data.frame(covariates)) {
          ## names of 'external' covariates to be found
          covnames <- variablesinformula(trend)
          if(want.subset)
            covnames <- union(covnames, all.vars(subsetexpr))
          if(allcovar)
            covnames <- union(covnames, names(covariates))
          covnames <- setdiff(covnames, c("x", "y", "marks"))
          ## resolve 'external' covariates
          tenv <- environment(trend)
          covariates <- getdataobjects(covnames, tenv, covariates, fatal=TRUE)
          updatecovariates <- any(attr(covariates, "external"))
        }
        ## extract values of covariates ('internal' and 'external')
        covariates.df <- mpl.get.covariates(covariates, P, Pname, covfunargs)
      }
      if(savecomputed)
        computed$covariates.df <- covariates.df
    } 

    ## Form the weights and the ``response variable''.

    if("dotmplbase" %in% names.precomputed) 
      .mpl <- precomputed$dotmplbase
    else {
      nQ <- n.quad(Q)
      wQ <- w.quad(Q)
      mQ <- marks.quad(Q)   ## is NULL for unmarked patterns
      zQ <- is.data(Q)
      yQ <- numeric(nQ)
      yQ[zQ] <- 1/wQ[zQ]
      zeroes <- attr(wQ, "zeroes")
      sQ <- if(is.null(zeroes)) rep.int(TRUE, nQ) else !zeroes
      ## tweak weights ONLY
      if(!is.null(weightfactor))
        wQ <- wQ * weightfactor
      ## pack up
      .mpl <- list(W      = wQ,
                   Z      = zQ,
                   Y      = yQ,
                   MARKS  = mQ, 
                   SUBSET = sQ)
    }

    if(savecomputed)
      computed$dotmplbase <- .mpl
  
    glmdata <- data.frame(.mpl.W = .mpl$W,
                          .mpl.Y = .mpl$Y)

    ## count data and dummy points in specified subset
    izdat <- .mpl$Z[.mpl$SUBSET]
    ndata <- sum(izdat)
#    ndummy <- sum(!izdat)
    
    ## Determine the domain of integration for the pseudolikelihood.
    if(correction == "border") {
      bdP <-
        if("bdP" %in% names.precomputed)
          precomputed$bdP
        else
          bdist.points(P)
      if(savecomputed)
        computed$bdP <- bdP
      .mpl$DOMAIN <- (bdP >= rbord)
    }

    skip.border <- skip.border && (correction == "border")
  
    ####################### T r e n d ##############################

    internal.names <- c(".mpl.W", ".mpl.Y", ".mpl.Z", ".mpl.SUBSET",
                        "SUBSET", ".mpl")

    reserved.names <- c("x", "y", "marks", internal.names)

    if(allcovar || want.trend || want.subset) {
      trendvariables <- variablesinformula(trend)
      ## Check for use of internal names in trend
      cc <- check.clashes(internal.names, trendvariables, "the model formula")
      if(cc != "") stop(cc)
      if(want.subset) {
        subsetvariables <- all.vars(subsetexpr)
        cc <- check.clashes(internal.names, trendvariables,
                            "the subset expression")
        if(cc != "") stop(cc)
        trendvariables <- union(trendvariables, subsetvariables)
      }
      ## Standard variables
      if(allcovar || "x" %in% trendvariables)
        glmdata <- data.frame(glmdata, x=P$x)
      if(allcovar || "y" %in% trendvariables)
        glmdata <- data.frame(glmdata, y=P$y)
      if(("marks" %in% trendvariables) || !is.null(.mpl$MARKS)) {
        if(is.null(.mpl$MARKS))
          stop("Model formula depends on marks, but data do not have marks",
               call.=FALSE)
        glmdata <- data.frame(glmdata, marks=.mpl$MARKS)
      }
      ##
      ## Check covariates
      if(!is.null(covariates.df)) {
        ## Check for duplication of reserved names
        cc <- check.clashes(reserved.names, names(covariates),
                            sQuote("covariates"))
        if(cc != "") stop(cc)
        ## Take only those covariates that are named in the trend formula
        if(!allcovar) 
          needed <- names(covariates.df) %in% trendvariables
        else
          needed <- rep.int(TRUE, ncol(covariates.df))
        if(any(needed)) {
          covariates.needed <- covariates.df[, needed, drop=FALSE]
          ##  Append to `glmdata'
          glmdata <- data.frame(glmdata,covariates.needed)
          ##  Ignore any quadrature points that have NA's in the covariates
          nbg <- is.na(covariates.needed)
          if(any(nbg)) {
            offending <- matcolany(nbg)
            covnames.na <- names(covariates.needed)[offending]
            quadpoints.na <- matrowany(nbg)
            n.na <- sum(quadpoints.na)
            n.tot <- length(quadpoints.na)
            errate <- n.na/n.tot
            pcerror <- round(signif(100 * errate, 2), 2)
            complaint <- paste("Values of the",
                               ngettext(length(covnames.na),
                                        "covariate", "covariates"),
                               paste(sQuote(covnames.na), collapse=", "),
                               "were NA or undefined at",
                               paste(pcerror, "%",
                                     " (", 
                                     n.na,
                                     " out of ",
                                     n.tot,
                                     ")",
                                     sep=""),
                               "of the", Pname)
            warning(paste(complaint,
                          ". Occurred while executing: ",
                          callstring, sep=""),
                    call. = FALSE)
            .mpl$SUBSET <-  .mpl$SUBSET & !quadpoints.na
            details <- list(covnames.na   = covnames.na,
                            quadpoints.na = quadpoints.na,
                            print         = complaint)
            problems <- append(problems,
                               list(na.covariates=details))
          }
        }
      }
    }

    ###################### I n t e r a c t i o n ####################

    Vnames <- NULL
    IsOffset <- NULL
    forbid <- NULL
    
    if(want.inter) {
      ## Form the matrix of "regression variables" V.
      ## The rows of V correspond to the rows of P (quadrature points)
      ## while the column(s) of V are the regression variables (log-potentials)

      E <- precomputed$E %orifnull% equalpairs.quad(Q)

      if(!skip.border) {
        ## usual case
        V <- evalInteraction(X, P, E, interaction, correction,
                             ...,
                             splitInf=splitInf,
                             precomputed=precomputed,
                             savecomputed=savecomputed)
      } else {
        ## evaluate only in eroded domain
        if(all(c("Esub", "Usub", "Retain") %in% names.precomputed)) {
          ## use precomputed data
          Psub <- precomputed$Usub
          Esub <- precomputed$Esub
          Retain <- precomputed$Retain
        } else {
          ## extract subset of quadrature points
          Retain <- .mpl$DOMAIN | is.data(Q)
          Psub <- P[Retain]
          ## map serial numbers in P to serial numbers in Psub
          Pmap <- cumsum(Retain)
          ## extract subset of equal-pairs matrix
          keepE <- Retain[ E[,2] ]
          Esub <- E[ keepE, , drop=FALSE]
          ## adjust indices in equal pairs matrix
          Esub[,2] <- Pmap[Esub[,2]]
        }
        ## call evaluator on reduced data
        if(all(c("X", "Q", "U") %in% names.precomputed)) {
          subcomputed <- resolve.defaults(list(E=Esub, U=Psub, Q=Q[Retain]),
                                          precomputed)
        } else subcomputed <- NULL
        if(clip.interaction) {
          ## normal
          V <- evalInteraction(X, Psub, Esub, interaction, correction,
                               ...,
                               splitInf=splitInf,
                               precomputed=subcomputed,
                               savecomputed=savecomputed)
        } else {
          ## ignore window when calculating interaction
          ## by setting 'W=NULL' (currently detected only by AreaInter)
          V <- evalInteraction(X, Psub, Esub, interaction, correction,
                               ...,
                               W=NULL,
                               splitInf=splitInf,
                               precomputed=subcomputed,
                               savecomputed=savecomputed)
        }
        if(savecomputed) {
          computed$Usub <- Psub
          computed$Esub <- Esub
          computed$Retain <- Retain
        }
      }
      
      if(!is.matrix(V))
        stop("interaction evaluator did not return a matrix")
      
      ## extract information about offsets
      IsOffset <- attr(V, "IsOffset")
      if(is.null(IsOffset)) IsOffset <- FALSE

      if(splitInf) {
        ## extract information about hard core terms
        forbid <- attr(V, "-Inf") %orifnull% logical(nrow(V))
      } 

      if(skip.border) {
        ## fill in the values in the border region with zeroes.
        Vnew <- matrix(0, nrow=npoints(P), ncol=ncol(V))
        colnames(Vnew) <- colnames(V)
        Vnew[Retain, ] <- V
        ## retain attributes
        attr(Vnew, "IsOffset") <- IsOffset
        attr(Vnew, "computed") <- attr(V, "computed")
        attr(Vnew, "POT") <- attr(V, "POT")
        V <- Vnew
        if(splitInf) {
          fnew <- logical(nrow(Vnew))
          fnew[Retain] <- forbid
          forbid <- fnew
        }
      }
    
      ## extract intermediate computation results 
      if(savecomputed)
        computed <- resolve.defaults(attr(V, "computed"), computed)

      ## Augment data frame by appending the regression variables
      ## for interactions.
      ##
      ## First determine the names of the variables
      ##
      Vnames <- dimnames(V)[[2]]
      if(is.null(Vnames)) {
        ## No names were provided for the columns of V.
        ## Give them default names.
        ## In ppm the names will be "Interaction"
        ##   or "Interact.1", "Interact.2", ...
        ## In mppm an alternative tag will be specified by vnamebase.
        nc <- ncol(V)
        Vnames <- if(nc == 1) vnamebase[1] else paste0(vnamebase[2], 1:nc)
        dimnames(V) <- list(dimnames(V)[[1]], Vnames)
      } else if(!is.null(vnameprefix)) {
        ## Variable names were provided by the evaluator (e.g. MultiStrauss).
        ## Prefix the variable names by a string
        ## (typically required by mppm)
        Vnames <- paste(vnameprefix, Vnames, sep="")
        dimnames(V) <- list(dimnames(V)[[1]], Vnames)
      }
      
      ## Check the names are valid as column names in a dataframe
      okVnames <- make.names(Vnames, unique=TRUE)
      if(any(Vnames != okVnames)) {
        warning(paste("Names of interaction terms",
                      "contained illegal characters;",
                      "names have been repaired."))
        Vnames <- okVnames
      }
    
      ##   Check for name clashes between the interaction variables
      ##   and the formula
      cc <- check.clashes(Vnames, termsinformula(trend), "model formula")
      if(cc != "") stop(cc)
      ##   and with the variables in 'covariates'
      if(!is.null(covariates)) {
        cc <- check.clashes(Vnames, names(covariates), sQuote("covariates"))
        if(cc != "") stop(cc)
      }

      ## OK. append variables.
      glmdata <- data.frame(glmdata, V)   

      ## check IsOffset matches Vnames
      if(length(IsOffset) != length(Vnames)) {
        if(length(IsOffset) == 1)
          IsOffset <- rep.int(IsOffset, length(Vnames))
        else
          stop("Internal error: IsOffset has wrong length", call.=FALSE)
      }
  
      ## Keep only those quadrature points for which the
      ## conditional intensity is nonzero. 

      ##KEEP  <- apply(V != -Inf, 1, all)
      .mpl$KEEP  <- matrowall(V != -Inf)

      .mpl$SUBSET <- .mpl$SUBSET & .mpl$KEEP

      ## Check that there are at least some data and dummy points remaining
      datremain <- .mpl$Z[.mpl$SUBSET]
      somedat <- any(datremain)
      somedum <- !all(datremain)
      if(warn.unidentifiable && !(somedat && somedum)) {
        ## Model would be unidentifiable if it were fitted.
        ## Register problem
        is.identifiable <- FALSE
        if(ndata == 0) {
          complaint <- "model is unidentifiable: data pattern is empty"
        } else {
          offending <- !c(somedat, somedum)
          offending <- c("all data points", "all dummy points")[offending]
          offending <- paste(offending, collapse=" and ")
          complaint <- paste("model is unidentifiable:",
                             offending, "have zero conditional intensity")
        }
        details <- list(data=!somedat,
                        dummy=!somedum,
                        print=complaint)
        problems <- append(problems, list(unidentifiable=details))
      }

      ## check whether the model has zero likelihood:
      ## check whether ANY data points have zero conditional intensity
      if(any(.mpl$Z & !.mpl$KEEP)) {
        howmany <- sum(.mpl$Z & !.mpl$KEEP)
        complaint <- paste(howmany,
                           "data point(s) are illegal",
                           "(zero conditional intensity under the model)")
        details <- list(illegal=howmany,
                        print=complaint)
        problems <- append(problems, list(zerolikelihood=details))
        if(warn.illegal && is.identifiable)
          warning(paste(complaint,
                        ". Occurred while executing: ",
                        callstring, sep=""),
                  call. = FALSE)
        likelihood.is.zero <- TRUE
      }
    }
  
    ##################     S u b s e t   ###################

    if(correction == "border") 
      .mpl$SUBSET <- .mpl$SUBSET & .mpl$DOMAIN
  
    if(!is.null(subsetexpr)) {
      ## user-defined subset expression
      USER.SUBSET <- eval(subsetexpr, glmdata, environment(trend))
      if(is.owin(USER.SUBSET)) {
        USER.SUBSET <- inside.owin(P$x, P$y, USER.SUBSET)
      } else if(is.im(USER.SUBSET)) {
        USER.SUBSET <- as.logical(USER.SUBSET[P, drop=FALSE])
        if(anyNA(USER.SUBSET)) 
          USER.SUBSET[is.na(USER.SUBSET)] <- FALSE
      }
      if(!(is.logical(USER.SUBSET) || is.numeric(USER.SUBSET)))
        stop("Argument 'subset' should yield logical values", call.=FALSE)
      if(anyNA(USER.SUBSET)) {
        USER.SUBSET[is.na(USER.SUBSET)] <- FALSE
        warning("NA values in argument 'subset' were changed to FALSE",
                call.=FALSE)
      }
      .mpl$SUBSET <- .mpl$SUBSET & USER.SUBSET
    }
                        
    glmdata <- cbind(glmdata,
                     data.frame(.mpl.SUBSET=.mpl$SUBSET,
                                stringsAsFactors=FALSE))

    #################  F o r m u l a   ##################################

    if(!want.trend) trend <- ~1 
    trendpart <- paste(as.character(trend), collapse=" ")
    if(!want.inter)
      rhs <- trendpart
    else {
      VN <- Vnames
      ## enclose offset potentials in 'offset(.)'
      if(any(IsOffset))
        VN[IsOffset] <- paste("offset(", VN[IsOffset], ")", sep="")
      rhs <- paste(c(trendpart, VN), collapse= "+")
    }
    fmla <- paste(".mpl.Y ", rhs)
    fmla <- as.formula(fmla)

    ##  character string of trend formula (without Vnames)
    trendfmla <- paste(".mpl.Y ", trendpart)

    ####
    result <- list(fmla=fmla, trendfmla=trendfmla,
                   covariates=if(updatecovariates) covariates else NULL,
                   glmdata=glmdata, Vnames=Vnames, IsOffset=IsOffset,
                   subsetexpr=subsetexpr,
                   problems=problems,
                   likelihood.is.zero=likelihood.is.zero,
                   is.identifiable=is.identifiable,
                   computed=computed,
                   vnamebase=vnamebase, vnameprefix=vnameprefix,
                   forbid=forbid)
    return(result)
  }

  check.clashes <- function(forbidden, offered, where) {
    name.match <- outer(forbidden, offered, "==")
    if(any(name.match)) {
      is.matched <- apply(name.match, 2, any)
      matched.names <- (offered)[is.matched]
      if(sum(is.matched) == 1) {
        return(paste("The variable",sQuote(matched.names),
                   "in", where, "is a reserved name"))
      } else {
        return(paste("The variables",
                   paste(sQuote(matched.names), collapse=", "),
                   "in", where, "are reserved names"))
      }
    }
    return("")
  }

  mpl.prepare
})



####################################################################
####################################################################

mpl.usable <- function(x) {
  ## silently remove covariates that don't have recognised format
  if(length(x) == 0 || is.data.frame(x)) return(x)
  isim   <- sapply(x, is.im)
  isfun  <- sapply(x, is.function)
  iswin  <- sapply(x, is.owin)
  istess <- sapply(x, is.tess)
  isnum  <- sapply(x, is.numeric) & (lengths(x) == 1)
  recognised <- isim | isfun | iswin | istess | isnum
  if(!all(recognised)) 
    x <- x[recognised]
  return(x)
}

mpl.get.covariates <- local({

  mpl.get.covariates <- function(covariates, locations, type="locations",
                                 covfunargs=list(),
                                 need.deriv=FALSE) {
    covargname <- sQuote(short.deparse(substitute(covariates)))
    locargname <- sQuote(short.deparse(substitute(locations)))
    if(is.null(covfunargs)) covfunargs <- list()
    ##
    x <- locations$x
    y <- locations$y
    if(is.null(x) || is.null(y)) {
      xy <- xy.coords(locations)
      x <- xy$x
      y <- xy$y
    }
    if(is.null(x) || is.null(y))
      stop(paste("Can't interpret", locargname, "as x,y coordinates"))
    n <- length(x)
    if(is.data.frame(covariates)) {
      if(nrow(covariates) != n)
        stop(paste("Number of rows in", covargname, 
                   "does not equal the number of", type))
      return(covariates)
    } else if(is.list(covariates)) {
      if(length(covariates) == 0)
        return(as.data.frame(matrix(, n, 0)))
      isim   <- unlist(lapply(covariates, is.im))
      isfun  <- unlist(lapply(covariates, is.function))
      iswin  <- unlist(lapply(covariates, is.owin))
      istess <- unlist(lapply(covariates, is.tess))
      isnum  <- unlist(lapply(covariates, is.number))
      if(!all(isim | isfun | isnum | iswin | istess))
        stop(paste("Each entry in the list", covargname, 
                   "should be an image, a function,",
                   "a window, a tessellation or a single number"))
      if(sum(nzchar(names(covariates))) < length(covariates))
        stop(paste("Some entries in the list",
                   covargname, "are un-named"))
      ## look up values of each covariate at the quadrature points
      values <- unclass(covariates)
      values[isim] <- lapply(covariates[isim], lookup.im, x=x, y=y,
                             naok=TRUE, strict=FALSE)
      values[isfun] <- vf <- lapply(covariates[isfun], evalfxy, x=x, y=y,
                                    extra=covfunargs)
      values[isnum] <- lapply(covariates[isnum], rep, length(x))
      values[iswin] <- lapply(covariates[iswin], insidexy, x=x, y=y)
      values[istess] <- lapply(covariates[istess], tileindex, x=x, y=y)
      result <- as.data.frame(values)
      if(need.deriv && any(isfun)) {
        ## check for gradient/hessian attributes of function values
        grad <- lapply(vf, attr, which="gradient")
        hess <- lapply(vf, attr, which="hessian")
        grad <- grad[!unlist(lapply(grad, is.null))]
        hess <- hess[!unlist(lapply(hess, is.null))]
        if(length(grad) > 0 || length(hess) > 0)
          attr(result, "derivatives") <- list(gradient=grad, hessian=hess)
      }
      return(result)
    } 
    stop(paste(covargname, "must be either a data frame or a list"))
  }

  ## functions for 'apply'
  evalfxy <- function(f, x, y, extra) {
    if(length(extra) == 0)
      return(f(x,y))
    ## extra arguments must be matched explicitly by name
    ok <- names(extra) %in% names(formals(f))
    z <- do.call(f, append(list(x,y), extra[ok]))
    return(z)
  }

  insidexy <- function(w, x, y) { inside.owin(x, y, w) }

  is.number <- function(x) { is.numeric(x) && (length(x) == 1) }

  mpl.get.covariates
})

bt.frame <- function(Q, trend=~1, interaction=NULL,
                      ...,
                      covariates=NULL,
                      correction="border", rbord=0,
                      use.gam=FALSE, allcovar=FALSE) {
  prep <- mpl.engine(Q, trend=trend, interaction=interaction,
                     ..., covariates=covariates,
                     correction=correction, rbord=rbord,
                     use.gam=use.gam, allcovar=allcovar,
                     preponly=TRUE, forcefit=TRUE)
  class(prep) <- c("bt.frame", class(prep))
  return(prep)
}


print.bt.frame <- function(x, ...) {
  cat("Model frame for Berman-Turner device\n")
  df <- x$glmdata
  cat(paste("$glmdata: Data frame with", nrow(df), "rows and",
            ncol(df), "columns\n"))
  cat("          Column names:\t")
  cat(paste(paste(names(df),collapse="\t"), "\n"))
  cat("Complete model formula ($fmla):\t")
  print(x$fmla)
  info <- x$info
  if(info$want.trend) {
    cat("Trend:\tyes\nTrend formula string ($trendfmla):\t")
    cat(paste(x$trendfmla, "\n"))
  } else cat("Trend:\tno\n")
  cat("Interaction ($info$interaction):\t")
  inte <- info$interaction
  if(is.null(inte))
    inte <- Poisson()
  print(inte, family=FALSE, brief=TRUE)
  if(!is.poisson.interact(inte)) {
    cat("Internal names of interaction variables ($Vnames):\t")
    cat(paste(x$Vnames, collapse="\t"))
    cat("\n")
  }
  edge <- info$correction
  cat(paste("Edge correction ($info$correction):\t", sQuote(edge), "\n"))
  if(edge == "border") 
    cat(paste("\tBorder width ($info$rbord):\t", info$rbord, "\n"))
  if(length(x$problems) > 0) {
    cat("Problems:\n")
    print(x$problems)
  }
  if(length(x$computed) > 0)
    cat(paste("Frame contains saved computations for",
              commasep(dQuote(names(x$computed)))))
  return(invisible(NULL))
}
  
partialModelMatrix <- function(X, D, model, callstring="", ...) {
  ## X = 'data'
  ## D = 'dummy'
  Q <- quad(X,D)
  P <- union.quad(Q)
  trend <- model$trend
  inter <- model$interaction
  covar <- model$covariates
  prep  <- mpl.prepare(Q, X, P, trend, inter, covar,
                       correction=model$correction,
                       rbord=model$rbord,
                       Pname="data points", callstring=callstring,
                       warn.unidentifiable=FALSE,
                       ...)
  fmla    <- prep$fmla
  glmdata <- prep$glmdata
  mof <- model.frame(fmla, glmdata)
  mom <- model.matrix(fmla, mof)

  modelnames <- names(coef(model))
  modelnames <- sub("log(lambda)", "(Intercept)", modelnames, fixed=TRUE)
  if(!isTRUE(all.equal(colnames(mom), modelnames)))
    warning(paste("Internal error: mismatch between",
                  "column names of model matrix",
                  "and names of coefficient vector in fitted model"))

  attr(mom, "mplsubset") <- glmdata$.mpl.SUBSET
  attr(mom, "-Inf") <- prep$forbid
  return(mom)
}
  
oversize.quad <- function(Q, ..., nU, nX, p=1) {
  ## Determine whether the quadrature scheme is
  ## too large to handle in one piece (in mpl)
  ## for a generic interaction
  ##    nU = number of quadrature points
  ##    nX = number of data points
  ##    p = dimension of statistic 
  if(missing(nU))
    nU <- n.quad(Q)
  if(missing(nX))
    nX <- npoints(Q$data)
  nmat <- as.double(nU) * nX
  nMAX <- spatstat.options("maxmatrix")/p
  needsplit <- (nmat > nMAX)
  return(needsplit)
}

quadBlockSizes <- function(nX, nD, p=1,
                           nMAX=spatstat.options("maxmatrix")/p,
                           announce=TRUE) {
  if(is.quad(nX) && missing(nD)) {
    nD <- npoints(nX$dummy)
    nX <- npoints(nX$data)
  }
  ## Calculate number of dummy points in largest permissible X * (X+D) matrix 
  nperblock <- max(1, floor(nMAX/nX - nX))
  ## determine number of such blocks 
  nblocks <- ceiling(nD/nperblock)
  ## make blocks roughly equal (except for the last one)
  nperblock <- min(nperblock, ceiling(nD/nblocks))
  ## announce
  if(announce && nblocks > 1) {
    msg <- paste("Large quadrature scheme",
                 "split into blocks to avoid memory size limits;",
                 nD, "dummy points split into",
                 nblocks, "blocks,")
    nfull <- nblocks - 1
    nlastblock <- nD - nperblock * nfull
    if(nlastblock == nperblock) {
      msg <- paste(msg,
                   "each containing",
                   nperblock, "dummy points")
    } else {
      msg <- paste(msg,
                   "the first",
                   ngettext(nfull, "block", paste(nfull, "blocks")),
                   "containing",
                   nperblock,
                   ngettext(nperblock, "dummy point", "dummy points"),
                   "and the last block containing",
                   nlastblock, 
                   ngettext(nlastblock, "dummy point", "dummy points"))
    }
    message(msg)
  } else nlastblock <- nperblock
  return(list(nblocks=nblocks, nperblock=nperblock, nlastblock=nlastblock))
}
    
## function that should be called to evaluate interaction terms
## between quadrature points and data points

evalInteraction <- function(X, P, E = equalpairs(P, X), 
                            interaction, correction,
                            splitInf=FALSE,
                            ...,
                            precomputed=NULL,
                            savecomputed=FALSE) {

  ## evaluate the interaction potential
  ## (does not assign/touch the variable names)

  verifyclass(interaction, "interact")

  ## handle Poisson case
  if(is.poisson(interaction)) {
    out <- matrix(numeric(0), nrow=npoints(P), ncol=0)
    attr(out, "IsOffset") <- logical(0)
    if(splitInf)
      attr(out, "-Inf") <- logical(nrow(out))
    return(out)
  }
  
  ## determine whether to use fast evaluation in C
  dofast <- (spatstat.options("fasteval") %in% c("on", "test")) &&
            !is.null(cando <- interaction$can.do.fast) &&
            cando(X, correction, interaction$par) &&
            !is.null(interaction$fasteval)

  ## determine whether to split quadscheme into blocks
  if(dofast) {
    dosplit <- FALSE
  } else {
    ## decide whether the quadrature scheme is too large to handle in one piece
    needsplit <- oversize.quad(nU=npoints(P), nX=npoints(X))

    ## not implemented when savecomputed=TRUE
    dosplit   <- needsplit && !savecomputed
    if(needsplit && savecomputed)
      warning(paste("Oversize quadscheme cannot be split into blocks",
                    "because savecomputed=TRUE;",
                    "memory allocation error may occur"))
  }
  
  if(!dosplit) {
    ## normal case
    V <- evalInterEngine(X=X, P=P, E=E,
                         interaction=interaction,
                         correction=correction,
                         splitInf=splitInf,
                         ...,
                         precomputed=precomputed,
                         savecomputed=savecomputed)
  } else {
    ## Too many quadrature points: split into blocks
    nX <- npoints(X)
    nP <- npoints(P)
    ## Determine which evaluation points are data points
    Pdata <- E[,2]
    ## hence which are dummy points
    Pall <- seq_len(nP)
    Pdummy <- if(length(Pdata) > 0) Pall[-Pdata] else Pall
    nD <- length(Pdummy)
    ## calculate block sizes
    bls <- quadBlockSizes(nX, nD, announce=TRUE)
    nblocks    <- bls$nblocks
    nperblock  <- bls$nperblock
    ##
    seqX <- seq_len(nX)
    EX <- cbind(seqX, seqX)
    ##
    for(iblock in 1:nblocks) {
      first <- min(nD, (iblock - 1) * nperblock + 1)
      last  <- min(nD, iblock * nperblock)
      ## extract dummy points  
      Di <- P[Pdummy[first:last]]
      Pi <- superimpose(X, Di, check=FALSE, W=X$window)
      ## evaluate potential
      Vi <- evalInterEngine(X=X, P=Pi, E=EX, 
                            interaction=interaction,
                            correction=correction,
                            splitInf=splitInf,
                            ...,
                            savecomputed=FALSE)
      Mi <- attr(Vi, "-Inf")
      if(iblock == 1) {
        V <- Vi
        M <- Mi
      } else {
        ## tack on the glm variables for the extra DUMMY points only
        V <- rbind(V, Vi[-seqX, , drop=FALSE])
        if(splitInf && !is.null(M))
          M <- c(M, Mi[-seqX])
      }
    }
    ## The first 'nX' rows of V contain values for X.
    ## The remaining rows of V contain values for dummy points.
    if(length(Pdata) == 0) {
      ## simply discard rows corresponding to data
      V <- V[-seqX, , drop=FALSE]
      if(splitInf && !is.null(M)) 
        M <- M[-seqX]
    } else {
      ## replace data in correct position
      ii <- integer(nP)
      ii[Pdata] <- seqX
      ii[Pdummy] <- (nX+1):nrow(V)
      V <- V[ii, , drop=FALSE]
      if(splitInf && !is.null(M))
        M <- M[ii]
    }
    attr(V, "-Inf") <- M
  }
  return(V)
}

## workhorse function that actually calls relevant code to evaluate interaction

evalInterEngine <- function(X, P, E, 
                            interaction, correction,
                            splitInf=FALSE, 
                            ...,
                            Reach = NULL,
                            precomputed=NULL,
                            savecomputed=FALSE) {

  ## fast evaluator (C code) may exist
  fasteval <- interaction$fasteval
  cando    <- interaction$can.do.fast
  par      <- interaction$par
  feopt    <- spatstat.options("fasteval")
  dofast   <- !is.null(fasteval) &&
              (is.null(cando) || cando(X, correction,par)) &&
              (feopt %in% c("on", "test")) &&
              (!splitInf || ("splitInf" %in% names(formals(fasteval))))
    
  V <- NULL
  if(dofast) {
    if(feopt == "test")
      message("Calling fasteval")
    V <- fasteval(X, P, E,
                  interaction$pot, interaction$par, correction,
                  splitInf=splitInf, ...)
  }
  if(is.null(V)) {
    ## use generic evaluator for family
    evaluate <- interaction$family$eval
    evalargs <- names(formals(evaluate))

    if(splitInf && !("splitInf" %in% evalargs))
      stop("Sorry, the", interaction$family$name, "interaction family",
           "does not support calculation of the positive part",
           call.=FALSE)

    if(is.null(Reach)) Reach <- reach(interaction)

    if("precomputed" %in% evalargs) {
      ## Use precomputed data
      ## version 1.9-3 onward (pairwise and pairsat families)
      V <- evaluate(X, P, E,
                    interaction$pot,
                    interaction$par,
                    correction=correction,
                    splitInf=splitInf,
                    ...,
                    Reach=Reach, 
                    precomputed=precomputed,
                    savecomputed=savecomputed)
    } else {
      ## Cannot use precomputed data
      ## Object created by earlier version of ppm
      ## or not pairwise/pairsat interaction
      V <- evaluate(X, P, E,
                    interaction$pot,
                    interaction$par,
                    correction=correction,
                    splitInf=splitInf,
                    ..., Reach=Reach)
    }
  }

  return(V)
}

deltasuffstat <- local({
  
  deltasuffstat <- function(model, ...,
                            restrict=c("pairs", "first", "none"),
			    dataonly=TRUE,
                            sparseOK=FALSE,
                            quadsub=NULL,
                            force=FALSE,
                            warn.forced=FALSE,
                            verbose=warn.forced,
                            use.special=TRUE) {
    stopifnot(is.ppm(model))
    sparsegiven <- !missing(sparseOK)
    restrict <- match.arg(restrict)
    
    if(dataonly) {
      X <- data.ppm(model)
      nX <- npoints(X)
    } else {
      X <- quad.ppm(model)
      if(!is.null(quadsub)) {
        z <- is.data(X)
        z[quadsub] <- FALSE
        if(any(z))
          stop("subset 'quadsub' must include all data points", call.=FALSE)
        X <- X[quadsub]
      }
      nX <- n.quad(X)
    }
    ncoef <- length(coef(model))
    inte <- as.interact(model)

    if(!sparseOK && exceedsMaxArraySize(nX, nX, ncoef)) {
      if(sparsegiven)
        stop("Array dimensions too large", call.=FALSE)
      warning("Switching to sparse array code", call.=FALSE)
      sparseOK <- TRUE
    }

    zeroes <- if(!sparseOK) array(0, dim=c(nX, nX, ncoef)) else
                            sparse3Darray(dims=c(nX, nX, ncoef))
    
    if(is.poisson(inte))
      return(zeroes)
    
    ## Get names of interaction terms in model (including offsets)
    f <- fitin(model)
    Inames <- f$Vnames
    IsOffset <- f$IsOffset
    hasInf <- !identical(inte$hasInf, FALSE)
    
    ## Offset terms do not contribute to sufficient statistic
    if(all(IsOffset) && !hasInf)
      return(zeroes)
    
    ## Nontrivial interaction terms must be computed.
    ## Look for member function $delta2 in the interaction
    v <- NULL
    v.is.full <- FALSE
    if(use.special) {
      ## Use specialised $delta2 for interaction,if available
      if(is.function(delta2 <- inte$delta2)) 
        v <- delta2(X, inte, model$correction, sparseOK=sparseOK)
      ## Use generic $delta2 for the family, if available
      if(is.null(v) && is.function(delta2 <- inte$family$delta2))
        v <- delta2(X, inte, model$correction, sparseOK=sparseOK)
    }
    ## no luck?
    if(is.null(v)) {
      if(!force)
        return(NULL)
      ## use brute force algorithm
      if(warn.forced)
        warning("Reverting to brute force to compute interaction terms",
                call.=FALSE)
      v <- if(dataonly) deltasufX(model, sparseOK, verbose=verbose) else
           deltasufQ(model, quadsub, sparseOK, verbose=verbose)
      v.is.full <- TRUE
    }
    
    ## extract hard core information
    deltaInf <- attr(v, "deltaInf")

    ## ensure 'v' is a 3D array
    if(length(dim(v)) != 3) {
      if(is.matrix(v)) {
        v <- array(v, dim=c(dim(v), 1))
      } else if(inherits(v, "sparseMatrix")) {
        v <- as.sparse3Darray(v)
      }
    }

    if(!sparseOK) {
      if(inherits(v, "sparse3Darray"))
        v <- as.array(v)
      if(inherits(deltaInf, "sparseMatrix"))
        deltaInf <- as.matrix(deltaInf)
    }

    if(restrict != "none") {
      ## kill contributions from points outside the domain of pseudolikelihood
      ## (e.g. points in the border region)
      use <- if(dataonly) getppmdatasubset(model) else
             if(is.null(quadsub)) getglmsubset(model) else
             getglmsubset(model)[quadsub]
      if(any(kill <- !use)) {
        switch(restrict,
 	       pairs = { v[kill,kill,] <- 0 },
	       first = { v[kill,,] <- 0 },
	       none = {})
        if(!is.null(deltaInf)) {
          switch(restrict,
                 pairs = { deltaInf[kill,kill] <- FALSE },
                 first = { deltaInf[kill,] <- FALSE },
                 none = {})
        }
      }
    }

    ## Make output array, with planes corresponding to model coefficients
    if(v.is.full) {
      ## Planes of 'v' already correspond to coefficients of model
      cnames <- names(coef(model))
      ## Remove any offset interaction terms
      ## (e.g. Hardcore interaction): these do not contribute to suff stat
      if(any(IsOffset)) {
        retain <- is.na(match(cnames, Inames[IsOffset]))
        v <- v[ , , retain, drop=FALSE]
        Inames <- Inames[!IsOffset]
      }
      result <- v
    } else {
      ## Planes of 'v' correspond to interaction terms only.
      ## Fill out the first order terms with zeroes
      result <- zeroes
      if(length(Inames) != dim(v)[3])
        stop(paste("Internal error: deltasuffstat:",
                   "number of planes of v =", dim(v)[3],
                   "!= number of interaction terms =", length(Inames)),
             call.=FALSE)
      ## Offset terms do not contribute to sufficient statistic
      if(any(IsOffset)) {
        v <- v[ , , !IsOffset, drop=FALSE]
        Inames <- Inames[!IsOffset]
      }
      ## Map planes of 'v' into coefficients
      Imap <- match(Inames, names(coef(model)))
      if(anyNA(Imap))
        stop(paste("Internal error: deltasuffstat:",
                   "cannot match interaction coefficients"))
      if(length(Imap) > 0) {
        ## insert 'v' into array
        result[ , , Imap] <- v
      }
    }
    ## pack up
    attr(result, "deltaInf") <- deltaInf
    return(result)
  }

  ## compute deltasuffstat using partialModelMatrix

  deltasufX <- function(model, sparseOK=FALSE, verbose=FALSE) {
    stopifnot(is.ppm(model))
    X <- data.ppm(model)
    hasInf <- !identical(model$interaction$hasInf, FALSE)
  
    nX <- npoints(X)
    p <- length(coef(model))

    m <- model.matrix(model, splitInf=hasInf)
    if(hasInf) {
      isInf <- attr(m, "-Inf")
      hasInf <- !is.null(isInf)
    }

    isdata <- is.data(quad.ppm(model))
    m <- m[isdata, ,drop=FALSE]
    if(hasInf)
      isInf <- isInf[isdata]
    ok <- getppmdatasubset(model)
    
    ## canonical statistic before and after deleting X[j]
    ## mbefore[ , i, j] = h(X[i] | X)
    ## mafter[ , i, j] = h(X[i] | X[-j])
    ## where h(u|x) is the canonical statistic of the *positive* cif
    dimwork <- c(p, nX, nX)
    if(!sparseOK) {
      mafter <- mbefore <- array(t(m), dim=dimwork)
      isInfafter <- isInfbefore <-
        if(!hasInf) NULL else matrix(isInf, dim=dimwork[-1]) 
    } else {
      ## make empty arrays; fill in values later
      ## (but only where they might change)
      mafter <- mbefore <- sparse3Darray(dims=dimwork)
      isInfafter <- isInfbefore <- if(!hasInf) NULL else
        sparseMatrix(i=integer(0), j=integer(0), x=logical(0),
                     dims=dimwork[-1])
    }

    ## identify close pairs
    R <- reach(model)
    if(is.finite(R)) {
      cl <- closepairs(X, R, what="indices")
      I <- cl$i
      J <- cl$j
      cl2 <- closepairs(X, 2*R, what="indices")
      I2 <- cl2$i
      J2 <- cl2$j
    } else {
      ## either infinite reach, or something wrong
      IJ <- expand.grid(I=1:nX, J=1:nX)
      IJ <- subset(IJ, I != J)
      I2 <- I <- IJ$I
      J2 <- J <- IJ$J
    }
    ## DO NOT RESTRICT - THIS IS NOW DONE IN deltasuffstat
    ## filter:  I and J must both belong to the nominated subset 
    ## okIJ <- ok[I] & ok[J]
    ## I <- I[okIJ]
    ## J <- J[okIJ]
    ##
    if(length(I) > 0 && length(J) > 0) {
      ## .............. loop over pairs ........................
      uniqueI <- unique(I)
      npairs <- length(uniqueI)
      pstate <- list()
      if(verbose) 
        splat("Examining", npairs, "pairs of data points...")
      ## The following ensures that 'empty' and 'X' have compatible marks 
      empty <- X[integer(0)]
      ## 
      ## Run through pairs
      for(iter in seq_len(npairs)) {
        i <- uniqueI[iter]
        ## all points within 2R
        J2i <- unique(J2[I2==i])
        ## all points within R
        Ji  <- unique(J[I==i])
        nJi <- length(Ji)
        if(nJi > 0) {
          Xi <- X[i]
          ## neighbours of X[i]
          XJi <- X[Ji]
          ## replace X[-i] by X[-i] \cap b(0, 2R)
          X.i <- X[J2i]
          nX.i <- length(J2i)
          ## index of XJi in X.i
          J.i <- match(Ji, J2i)
          if(anyNA(J.i))
            stop("Internal error: Ji not a subset of J2i")
          ## values of sufficient statistic 
          ##    h(X[j] | X[-i]) = h(X[j] | X[-c(i,j)]
          ## for all j
          pmj <- snipModelMatrix(X.i, empty, model, J.i, hasInf)
          if(hasInf) zzj <- attr(pmj, "-Inf")
          ## sufficient statistic in reverse order
          ##    h(X[i] | X[-j]) = h(X[i] | X[-c(i,j)]
          ## for all j
          pmi <- matrix(, nJi, p)
          zzi <- logical(nJi)
          for(k in 1:nJi) {
            ## j <- Ji[k]
            ## X.ij <- X[-c(i,j)]
            X.ij <- X.i[-J.i[k]]
            pmik <- snipModelMatrix(X.ij, Xi, model, nX.i, hasInf)
            pmi[k, ] <- pmik
            if(hasInf) zzi[k] <- attr(pmik, "-Inf")
          }
          ##
          if(!sparseOK) {
            mafter[ , Ji, i] <- t(pmj)
            mafter[ , i, Ji] <- t(pmi)
            if(hasInf) {
              isInfafter[Ji, i] <- zzj
              isInfafter[i, Ji] <- zzi
            }
          } else {
            mafter[ , Ji, i] <- array(t(pmj), dim=c(p, nJi, 1))
            mafter[ , i, Ji] <- array(t(pmi), dim=c(p, 1, nJi))
            mbefore[ , Ji, i] <- array(t(m[Ji,]), dim=c(p, nJi, 1))
            mbefore[ , i, Ji] <- array(m[i,], dim=c(p, 1, nJi))
            if(hasInf) {
              isInfafter[Ji, i] <- zzj
              isInfafter[i, Ji] <- zzi
              isInfbefore[Ji, i] <- isInf[Ji]
              isInfbefore[i, Ji] <- isInf[i]
            }
          } 
        }
        if(verbose)
          pstate <- progressreport(iter, npairs, state=pstate)
      }
    }
        
    ##  delta[ ,i,j] = h(X[i] | X) - h(X[i] | X[-j])
    delta <- mbefore - mafter
    ## delta[i, j, ] = h(X[i] | X) - h(X[i] | X[-j])
    delta <- aperm(delta, c(2,3,1))
    ##
    if(hasInf) {
      deltaInf <- isInfbefore - isInfafter
      attr(delta, "deltaInf") <- deltaInf
    }
    return(delta)
  }

  deltasufQ <- function(model, quadsub, sparseOK, verbose=FALSE) {
    stopifnot(is.ppm(model))
    hasInf <- !identical(model$interaction$hasInf, FALSE)

    p <- length(coef(model))
    
    Q <- quad.ppm(model)
    ok <- getglmsubset(model)

    m <- model.matrix(model, splitInf=hasInf)
    if(hasInf) {
      isInf <- attr(m, "-Inf")
      hasInf <- !is.null(isInf)
    }

    if(!is.null(quadsub)) {
      Q <- Q[quadsub]
      m <- m[quadsub, , drop=FALSE]
      ok <- ok[quadsub]
      if(hasInf) isInf <- isInf[quadsub]
    }
    
    X <- Q$data
    U <- union.quad(Q)
    nU <- npoints(U)
    nX <- npoints(X)
    isdata <- is.data(Q)
    isdummy <- !isdata
    
    ## canonical statistic before and after adding/deleting U[j]
    dimwork <- c(p, nU, nU)
    if(!sparseOK) {
      mafter <- mbefore <- array(t(m), dim=dimwork)
      delta <- array(0, dim=dimwork)
      isInfafter <- isInfbefore <- deltaInf <-
        if(!hasInf) NULL else matrix(isInf, dim=dimwork[-1])
    } else {
      ## make empty arrays; fill in values later
      ## [but only where they might change]
      mafter <- mbefore <- delta <- sparse3Darray(dims=dimwork)
      isInfafter <- isInfbefore <- deltaInf <-
        if(!hasInf) NULL else sparseMatrix(i=integer(0), j=integer(0),
                                           x=logical(0),
                                           dims=dimwork[-1])
    }
    ##   mbefore[ , i, j] = h(U[i] | X)
    ## For data points X[j]
    ##   mafter[ , i, j] = h(U[i] | X[-j])
    ##   delta[ , i, j] = h(U[i] | X) - h(U[i] | X[-j])
    ## For dummy points X[j]
    ##   mafter[ , i, j] = h(U[i] | X \cup U[j])
    ##   delta[ , i, j] = h(U[i] | X \cup U[j]) - h(U[i] | X)

    changesign <- ifelseAB(isdata, -1, 1)
  
    ## identify close pairs of quadrature points
    R <- reach(model)
    if(is.finite(R)) {
      cl <- closepairs(U, R, what="indices")
      I <- cl$i
      J <- cl$j
      cl2 <- closepairs(U, 2*R, what="indices")
      I2 <- cl2$i
      J2 <- cl2$j
    } else {
      ## either infinite reach, or something wrong
      IJ <- expand.grid(I=1:nU, J=1:nX)
      IJ <- IJ[ with(IJ, I != J), ]
      I2 <- I <- IJ$I
      J2 <- J <- IJ$J
    }

    ## filter:  I and J must both belong to the nominated subset 
    okIJ <- ok[I] & ok[J]
    I <- I[okIJ]
    J <- J[okIJ]
    ##
    if(length(I) > 0 && length(J) > 0) {
      ## .............. loop over pairs of quadrature points ...............
      ## Run through pairs
      uI <- unique(I)
      zI <- isdata[uI]
      uIdata <- uI[zI]
      uIdummy <- uI[!zI]
      nuIdata <- length(uIdata)
      nuIdummy <- length(uIdummy)
      if(verbose) 
        splat("Examining", nuIdata, "+", nuIdummy, "=", nuIdata + nuIdummy,
              "pairs of points")
      ## Run through pairs i, j where 'i' is a data point
      pstate <- list()
      for(iter in seq_len(nuIdata)) {
        i <- uIdata[iter]
        ## all DATA points within 2R of X[i]
        ## This represents X[-i] 
        J2i <- unique(J2[I2==i])
        J2i <- J2i[isdata[J2i]]
        ## all QUADRATURE points within R of X[i]
        Ji  <- unique(J[I==i])
        nJi <- length(Ji)
        if(nJi > 0) {
          isd <- isdata[Ji]
          ## data points which are neighbours of X[i]
          XJi <- X[Ji[isd]]
          ## dummy points which are neighbours of X[i]
          DJi <- U[Ji[!isd]]
          ## replace X[-i] by X[-i] \cap b(0, 2R)
          X.i <- X[J2i]
          nX.i <- length(J2i)
          ## index of XJi in X.i 
          J.i <- match(Ji[isd], J2i)
          if(anyNA(J.i))
            stop("Internal error: Ji[isd] not a subset of J2i")
          ## index of DJi in superimpose(X.i, DJi)
          JDi <- nX.i + seq_len(sum(!isd))
          ## values of sufficient statistic 
          ##    h(X[j] | X[-i]) = h(X[j] | X[-c(i,j)]
          ## for all j
          pmj <- snipModelMatrix(X.i, DJi, model, c(J.i, JDi), hasInf)
          ##
          mafter[ , Ji, i] <- t(pmj)
          if(hasInf)
            isInfafter[Ji, i] <- attr(pmj, "-Inf")
          if(sparseOK) {
            mbefore[ , Ji, i] <- array(t(m[Ji,]), dim=c(p, nJi, 1))
            if(hasInf) isInfbefore[Ji, i] <- isInf[Ji]
          }
        }
        if(verbose)
          pstate <- progressreport(iter, nuIdata, state=pstate)
      }
      ## Run through pairs i, j where 'i' is a dummy point
      pstate <- list()
      for(iter in seq_len(nuIdummy)) {
        i <- uIdummy[iter]
        ## all DATA points within 2R of U[i]
        J2i <- unique(J2[I2==i])
        J2i <- J2i[isdata[J2i]]
        ## all QUADRATURE points within R of U[i]
        Ji  <- unique(J[I==i])
        nJi <- length(Ji)
        if(nJi > 0) {
          isd <- isdata[Ji]
          JiData <- Ji[isd]
          JiDummy <- Ji[!isd]
          ## data points which are neighbours of U[i]
          XJi <- X[JiData]
          ## dummy points which are neighbours of U[i]
          DJi <- U[JiDummy]
          ## replace X \cup U[i] by (X \cap b(0, 2R)) \cup U[i]
          J2Ui <- c(J2i, i)
          XUi <- U[J2Ui]
          nXUi <- length(J2Ui)
          ## index of XJi in X.i 
          J.i <- match(JiData, J2Ui)
          if(anyNA(J.i))
            stop("Internal error: Ji[isd] not a subset of J2i")
          ## index of DJi in superimpose(X.i, DJi)
          JDi <- nXUi + seq_len(length(JiDummy))
          ## values of sufficient statistic 
          ##    h(X[j] | X \cup U[i]) 
          ## for all j
          pmj <- snipModelMatrix(XUi, DJi, model, c(J.i, JDi), hasInf)
          ##
          JiSort <- c(JiData, JiDummy)
          if(!sparseOK) {
            mafter[ , JiSort, i] <- t(pmj)
            if(hasInf)
              isInfafter[JiSort, i] <- attr(pmj, "-Inf")
          } else {
            mafter[ , JiSort, i]  <- array(t(pmj), dim=c(p, nJi, 1))
            mbefore[ , JiSort, i] <- array(t(m[JiSort,]), dim=c(p, nJi, 1))
            if(hasInf) {
              isInfafter[JiSort, i] <- attr(pmj, "-Inf")
              isInfbefore[JiSort, i] <- isInf[JiSort]
            }
          }
        }
        if(verbose)
          pstate <- progressreport(iter, nuIdummy, state=pstate)
      }
    }
        
    ##  delta[ ,i,j] = h(X[i] | X) - h(X[i] | X[-j])
    delta[ , , isdata] <- mbefore[, , isdata] - mafter[ , , isdata]
    ##  delta[ ,i,j] = h(X[i] | X \cup U[j]) - h(X[i] | X)
    delta[ , , isdummy] <- mafter[, , isdummy] - mbefore[ , , isdummy]
    ## rearrange: new delta[i,j,] = old delta[, i, j]
    delta <- aperm(delta, c(2,3,1))
    ##
    if(hasInf) {
      deltaInf[ , isdata]  <- isInfbefore[ , isdata] - isInfafter[ , isdata]
      deltaInf[ , isdummy] <- isInfafter[ , isdummy] - isInfbefore[ , isdummy]
      attr(delta, "deltaInf") <- deltaInf
    }
    return(delta)
  }

  snipModelMatrix <- function(X, D, model, retain, splitInf=FALSE) {
    M <- partialModelMatrix(X, D, model, splitInf=splitInf)
    if(splitInf) isInf <- attr(M, "-Inf")
    M <- M[retain, , drop=FALSE]
    if(splitInf) attr(M, "-Inf") <- isInf[retain]
    return(M)
  }
      
  deltasuffstat
})

