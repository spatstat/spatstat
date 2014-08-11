#    mpl.R
#
#	$Revision: 5.154 $	$Date: 2013/04/25 06:37:43 $
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

"mpl" <- function(Q,
         trend = ~1,
	 interaction = NULL,
         data = NULL,
	 correction="border",
	 rbord = 0,
         use.gam=FALSE) {
   .Deprecated("ppm", package="spatstat")
   ppm(Q=Q, trend=trend, interaction=interaction,
       covariates=data, correction=correction, rbord=rbord,
       use.gam=use.gam, method="mpl")
}

"mpl.engine" <- 
function(Q,
         trend = ~1,
	 interaction = NULL,
         ...,
         covariates = NULL,
         covfunargs = list(),
	 correction="border",
	 rbord = 0,
         use.gam=FALSE,
         gcontrol=list(),
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
         justQ = FALSE
) {
#
# Extract precomputed data if available
#  
  if(!is.null(precomputed$Q)) {
    Q <- precomputed$Q
    X <- precomputed$X
    P <- precomputed$U
  } else {
#
# Determine quadrature scheme from argument Q
#
    if(verifyclass(Q, "quad", fatal=FALSE)) {
      # user-supplied quadrature scheme - validate it
      validate.quad(Q, fatal=TRUE, repair=FALSE, announce=TRUE)
      # Extract data points
      X <- Q$data
    } else if(verifyclass(Q, "ppp", fatal = FALSE)) {
      # point pattern - create default quadrature scheme
      X <- Q
      Q <- quadscheme(X, nd=nd, eps=eps)
    } else 
      stop("First argument Q should be a point pattern or a quadrature scheme")
    
#
#
# Data and dummy points together
    P <- union.quad(Q)

  }
#
# secret exit  
  if(justQ) return(Q)
#
#  
  computed <- if(savecomputed) list(X=X, Q=Q, U=P) else NULL
#
# Validate main arguments
  if(!is.null(trend) && !inherits(trend, "formula"))
    stop(paste("Argument", sQuote("trend"), "must be a formula"))
  if(!is.null(interaction) && !inherits(interaction, "interact"))
    stop(paste("Argument", sQuote("interaction"), "has incorrect format"))
#
  check.1.real(rbord, "In ppm")
  explain.ifnot(rbord >= 0, "In ppm")
# rbord applies only to border correction
  if(correction != "border") rbord <- 0 
#
# Self-starting interaction?
#
  if(!is.null(ss <- interaction$selfstart)) 
    interaction <- ss(X, interaction)
  
#
# Interpret the call
want.trend <- !is.null(trend) && !identical.formulae(trend, ~1)
want.inter <- !is.null(interaction) && !is.null(interaction$family)
trend.formula <- if(want.trend) trend else (~1)

  
# Stamp with spatstat version number
  
spv <- package_version(versionstring.spatstat())
the.version <- list(major=spv$major,
                    minor=spv$minor,
                    release=spv$patchlevel,
                    date="$Date: 2013/04/25 06:37:43 $")

if(want.inter) {
  # ensure we're using the latest version of the interaction object
  if(outdated.interact(interaction)) 
    interaction <- update(interaction)
}

#  
  
if(!want.trend && !want.inter && !forcefit && !allcovar) {
  # the model is the uniform Poisson process
  # The MPLE (= MLE) can be evaluated directly
  npts <- npoints(X)
  W    <- as.owin(X)
  if(correction == "border" && rbord > 0) {
    W <- erosion(W, rbord)
    npts <- npoints(X[W])
  }
  volume <- area.owin(W) * markspace.integral(X)
  lambda <- npts/volume
  # fitted canonical coefficient
  co <- log(lambda)
  # asymptotic variance of canonical coefficient
  varcov <- matrix(1/npts, 1, 1)
  fisher <- matrix(npts,   1, 1)
  se <- sqrt(1/npts)
  # give names
  tag <- if(rename.intercept) "log(lambda)" else "(Intercept)"
  names(co) <- tag
  dimnames(varcov) <- dimnames(fisher) <- list(tag, tag)
  # maximised log likelihood
  maxlogpl <- if(npts == 0) 0 else npts * (log(lambda) - 1)
  #
  rslt <- list(
               method      = "mpl",
               fitter      = "exact",
               projected   = FALSE,
               coef        = co,
               trend       = NULL,
               interaction = NULL,
               fitin       = fii(),
               Q           = Q,
               maxlogpl    = maxlogpl,
               internal    = list(computed=computed, se=se),
               covariates  = covariates,  # covariates are still retained!
               covfunargs  = covfunargs,
	       correction  = correction,
               rbord       = rbord,
               terms       = terms(trend.formula),
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
                      allcovar=allcovar,
                      precomputed=precomputed, savecomputed=savecomputed,
                      covfunargs=covfunargs,
                      ...)
  # back door
if(preponly) {
  # exit now, returning prepared data frame and internal information
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
computed <- append(computed, prep$computed)
IsOffset <- prep$IsOffset

  
################# F i t    i t   ####################################

if(!is.identifiable) 
  stop(paste("in", callstring, ":", problems$unidentifiable$print),
       call.=FALSE)
  
# to avoid problem with package checker  
.mpl.W <- glmdata$.mpl.W
.mpl.SUBSET <- glmdata$.mpl.SUBSET

# determine algorithm control parameters
  if(is.null(gcontrol)) gcontrol <- list() else stopifnot(is.list(gcontrol))
  gc <- if(use.gam) "gam.control" else "glm.control"
  gcontrol <- do.call(gc, gcontrol)
  
# Fit the generalized linear/additive model.

if(is.null(famille)) {
  # the sanctioned technique, using `quasi' family
  if(want.trend && use.gam)
    FIT  <- gam(fmla, family=quasi(link=log, variance=mu), weights=.mpl.W,
                data=glmdata, subset=.mpl.SUBSET,
                control=gcontrol)
  else
    FIT  <- glm(fmla, family=quasi(link=log, variance=mu), weights=.mpl.W,
                data=glmdata, subset=.mpl.SUBSET,
                control=gcontrol, model=FALSE)
} else {
  # for experimentation only!
  if(is.function(famille))
    famille <- famille()
  stopifnot(inherits(famille, "family"))
  if(want.trend && use.gam)
    FIT  <- gam(fmla, family=famille, weights=.mpl.W,
                data=glmdata, subset=.mpl.SUBSET,
                control=gcontrol)
  else
    FIT  <- glm(fmla, family=famille, weights=.mpl.W,
                data=glmdata, subset=.mpl.SUBSET,
                control=gcontrol, model=FALSE)
}
  environment(FIT$terms) <- sys.frame(sys.nframe())

  
################  I n t e r p r e t    f i t   #######################

# Fitted coefficients

  co <- FIT$coef

# glm covariates
  W <- glmdata$.mpl.W
  SUBSET <- glmdata$.mpl.SUBSET        
  Z <- is.data(Q)
  Vnames <- prep$Vnames
        
# attained value of max log pseudolikelihood
  maxlogpl <-
    if(likelihood.is.zero) { -Inf } else 
    -(deviance(FIT)/2 + sum(log(W[Z & SUBSET])) + sum(Z & SUBSET))

# fitted interaction object
  fitin <- if(want.inter) fii(interaction, co, Vnames, IsOffset) else fii()
######################################################################
# Clean up & return 

rslt <- list(
             method       = "mpl",
             fitter       = if(use.gam) "gam" else "glm",
             projected    = FALSE,
             coef         = co,
             trend        = if(want.trend) trend       else NULL,
             interaction  = if(want.inter) interaction else NULL,
             fitin        = fitin,
             Q            = Q,
             maxlogpl     = maxlogpl, 
             internal     = list(glmfit=FIT, glmdata=glmdata, Vnames=Vnames,
                              IsOffset=IsOffset, fmla=fmla, computed=computed),
             covariates   = covariates,
             covfunargs   = covfunargs,
             correction   = correction,
             rbord        = rbord,
             terms        = terms(trend.formula),
             version      = the.version,
             problems     = problems)
class(rslt) <- "ppm"
return(rslt)
}  


##########################################################################
### /////////////////////////////////////////////////////////////////////
##########################################################################

mpl.prepare <- function(Q, X, P, trend, interaction, covariates, 
                        want.trend, want.inter, correction, rbord,
                        Pname="quadrature points", callstring="",
                        ...,
                        covfunargs=list(),
                        allcovar=FALSE,
                        precomputed=NULL, savecomputed=FALSE,
                        vnamebase=c("Interaction", "Interact."),
                        vnameprefix=NULL,
                        warn.illegal=TRUE,
                        warn.unidentifiable=TRUE) {
# Q: quadrature scheme
# X = data.quad(Q)
# P = union.quad(Q)
  
  if(missing(want.trend))
    want.trend <- !is.null(trend) && !identical.formulae(trend, ~1)
  if(missing(want.inter))
    want.inter <- !is.null(interaction) && !is.null(interaction$family)
    
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

# Extract covariate values
  if((allcovar || want.trend) && !is.null(covariates)) {
    if("covariates.df" %in% names.precomputed)
      covariates.df <- precomputed$covariates.df
    else 
      covariates.df <- mpl.get.covariates(covariates, P, Pname, covfunargs)
    if(savecomputed)
      computed$covariates.df <- covariates.df
  }
        
### Form the weights and the ``response variable''.

  if("dotmplbase" %in% names.precomputed) 
    .mpl <- precomputed$dotmplbase
  else {
    .mpl <- list()
    .mpl$W <- w.quad(Q)
    .mpl$Z <- is.data(Q)
    .mpl$Y <- .mpl$Z/.mpl$W
    .mpl$MARKS <- marks.quad(Q)  # is NULL for unmarked patterns
    n <- n.quad(Q)
    .mpl$SUBSET <- rep.int(TRUE, n)
	
    zeroes <- attr(.mpl$W, "zeroes")
    if(!is.null(zeroes))
      .mpl$SUBSET <-  !zeroes
  }

  if(savecomputed)
    computed$dotmplbase <- .mpl
  
  glmdata <- data.frame(.mpl.W = .mpl$W,
                        .mpl.Y = .mpl$Y)

  # count data and dummy points in specified subset
  izdat <- .mpl$Z[.mpl$SUBSET]
  ndata <- sum(izdat)
  ndummy <- sum(!izdat)
    
####################### T r e n d ##############################

  internal.names <- c(".mpl.W", ".mpl.Y", ".mpl.Z", ".mpl.SUBSET",
                        "SUBSET", ".mpl")

  reserved.names <- c("x", "y", "marks", internal.names)

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
  
  if(allcovar || want.trend) {
    trendvariables <- variablesinformula(trend)
    # Check for use of internal names in trend
    cc <- check.clashes(internal.names, trendvariables, "the model formula")
    if(cc != "") stop(cc)
    # Standard variables
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
    #
    # Check covariates
    if(!is.null(covariates)) {
    # Check for duplication of reserved names
    cc <- check.clashes(reserved.names, names(covariates),
                        sQuote("covariates"))
    if(cc != "") stop(cc)
    # Take only those covariates that are named in the trend formula
    if(!allcovar) 
      needed <- names(covariates.df) %in% trendvariables
    else
      needed <- rep.int(TRUE, ncol(covariates.df))
    if(any(needed)) {
      covariates.needed <- covariates.df[, needed, drop=FALSE]
    #  Append to `glmdata'
      glmdata <- data.frame(glmdata,covariates.needed)
    #  Ignore any quadrature points that have NA's in the covariates
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

  if(want.inter) {

    # Form the matrix of "regression variables" V.
    # The rows of V correspond to the rows of P (quadrature points)
    # while the column(s) of V are the regression variables (log-potentials)

    E <- equalpairs.quad(Q)

    V <- evalInteraction(X, P, E, interaction, correction,
                         ...,
                         precomputed=precomputed,
                         savecomputed=savecomputed)
    if(!is.matrix(V))
      stop("interaction evaluator did not return a matrix")
    # extract intermediate computation results 
    if(savecomputed)
      computed <- append(computed, attr(V, "computed"))

    # extract information about offsets
    IsOffset <- attr(V, "IsOffset")
    if(is.null(IsOffset)) IsOffset <- FALSE 
  
    # Augment data frame by appending the regression variables for interactions.
    #
    # First determine the names of the variables
    #
    Vnames <- dimnames(V)[[2]]
    if(is.null(Vnames)) {
      # No names were provided for the columns of V.
      # Give them default names.
      # In ppm the names will be "Interaction" or "Interact.1", "Interact.2", ...
      # In mppm an alternative tag will be specified by vnamebase.
      nc <- ncol(V)
      Vnames <- if(nc == 1) vnamebase[1] else paste(vnamebase[2], 1:nc, sep="")
      dimnames(V) <- list(dimnames(V)[[1]], Vnames)
    } else if(!is.null(vnameprefix)) {
      # Variable names were provided by the evaluator (e.g. MultiStrauss).
      # Prefix the variable names by a string
      # (typically required by mppm)
      Vnames <- paste(vnameprefix, Vnames, sep="")
      dimnames(V) <- list(dimnames(V)[[1]], Vnames)
    }

    # Check the names are valid as column names in a dataframe
    okVnames <- make.names(Vnames, unique=TRUE)
    if(any(Vnames != okVnames)) {
      warning("Names of interaction terms contained illegal characters; names have been repaired.")
      Vnames <- okVnames
    }
    
    #   Check for name clashes between the interaction variables
    #   and the formula
    cc <- check.clashes(Vnames, termsinformula(trend), "model formula")
    if(cc != "") stop(cc)
    #   and with the variables in 'covariates'
    if(!is.null(covariates)) {
      cc <- check.clashes(Vnames, names(covariates), sQuote("covariates"))
      if(cc != "") stop(cc)
    }

    # OK. append variables.
    glmdata <- data.frame(glmdata, V)   

    # check IsOffset matches Vnames
    if(length(IsOffset) != length(Vnames)) {
      if(length(IsOffset) == 1)
        IsOffset <- rep.int(IsOffset, length(Vnames))
      else
        stop("Internal error: IsOffset has wrong length", call.=FALSE)
    }
  
    # Keep only those quadrature points for which the
    # conditional intensity is nonzero. 

    #KEEP  <- apply(V != -Inf, 1, all)
    .mpl$KEEP  <- matrowall(V != -Inf)

    .mpl$SUBSET <- .mpl$SUBSET & .mpl$KEEP

    # Check that there are at least some data and dummy points remaining
    datremain <- .mpl$Z[.mpl$SUBSET]
    somedat <- any(datremain)
    somedum <- !all(datremain)
    if(warn.unidentifiable && !(somedat && somedum)) {
      # Model would be unidentifiable if it were fitted.
      # Register problem
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

    # check whether the model has zero likelihood:
    # check whether ANY data points have zero conditional intensity
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
##################   D a t a    f r a m e   ###################

# Determine the domain of integration for the pseudolikelihood.

  if(correction == "border") {

    if("bdP" %in% names.precomputed)
      bd <- precomputed$bdP
    else
      bd <- bdist.points(P)

    if(savecomputed)
      computed$bdP <- bd
  
    .mpl$DOMAIN <- (bd >= rbord)
    .mpl$SUBSET <- .mpl$DOMAIN & .mpl$SUBSET
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
    # enclose offset potentials in 'offset(.)'
    if(any(IsOffset))
      VN[IsOffset] <- paste("offset(", VN[IsOffset], ")", sep="")
    rhs <- paste(c(trendpart, VN), collapse= "+")
  }
  fmla <- paste(".mpl.Y ", rhs)
  fmla <- as.formula(fmla)

  ####  character string of trend formula (without Vnames)
  trendfmla <- paste(".mpl.Y ", trendpart)

#### 

return(list(fmla=fmla, trendfmla=trendfmla,
            glmdata=glmdata, Vnames=Vnames, IsOffset=IsOffset,
            problems=problems,
            likelihood.is.zero=likelihood.is.zero,
            is.identifiable=is.identifiable,
            computed=computed))

}



####################################################################
####################################################################

mpl.get.covariates <- function(covariates, locations, type="locations",
                               covfunargs=list()) {
  covargname <- sQuote(short.deparse(substitute(covariates)))
  locargname <- sQuote(short.deparse(substitute(locations)))
  if(is.null(covfunargs)) covfunargs <- list()
  #
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
    is.number <- function(x) { is.numeric(x) && (length(x) == 1) }
    isim  <- unlist(lapply(covariates, is.im))
    isfun <- unlist(lapply(covariates, is.function))
    iswin <- unlist(lapply(covariates, is.owin))
    isnum <- unlist(lapply(covariates, is.number))
    if(!all(isim | isfun | isnum | iswin))
      stop(paste("Each entry in the list", covargname, 
                 "should be an image, a function, a window or a single number"))
    if(sum(nzchar(names(covariates))) < length(covariates))
      stop(paste("Some entries in the list",
                 covargname, "are un-named"))
    # look up values of each covariate at the quadrature points
    values <- covariates
    evalfxy <- function(f, x, y, extra) {
      if(length(extra) == 0)
        return(f(x,y))
      # extra arguments must be matched explicitly by name
      ok <- names(extra) %in% names(formals(f))
      z <- do.call(f, append(list(x,y), extra[ok]))
      return(z)
    }
    insidexy <- function(w, x, y) { inside.owin(x, y, w) }
    values[isim] <- lapply(covariates[isim], lookup.im, x=x, y=y,
                           naok=TRUE, strict=FALSE)
    values[isfun] <- lapply(covariates[isfun], evalfxy, x=x, y=y,
                            extra=covfunargs)
    values[isnum] <- lapply(covariates[isnum], rep, length(x))
    values[iswin] <- lapply(covariates[iswin], insidexy, x=x, y=y)
    return(as.data.frame(values))
  } else
    stop(paste(covargname, "must be either a data frame or a list"))
}

bt.frame <- function(Q, trend=~1, interaction=NULL,
                      ...,
                      covariates=NULL,
                      correction="border", rbord=0,
                      use.gam=FALSE, allcovar=FALSE) {
  prep <- mpl.engine(Q=Q, trend=trend, interaction=interaction,
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
  if(edge == "border") {
    if(info$rbord==0)
      edge <- "none"
    else
      edge <- paste("border", paren(paste("rbord=", info$rbord)))
  }
  cat(paste("Edge correction ($info$correction):\t", edge, "\n"))
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
  # X = 'data'
  # D = 'dummy'
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

  if(!identical(all.equal(colnames(mom), names(coef(model))), TRUE))
    warning("Internal error: mismatch between column names of model matrix and names of coefficient vector in fitted model")

  attr(mom, "mplsubset") <- glmdata$.mpl.SUBSET
  return(mom)
}
  
oversize.quad <- function(Q, ..., nU, nX) {
  # Determine whether the quadrature scheme is
  # too large to handle in one piece (in mpl)
  # for a generic interaction
  #    nU = number of quadrature points
  #    nX = number of data points
  if(missing(nU))
    nU <- n.quad(Q)
  if(missing(nX))
    nX <- npoints(Q$data)
  nmat <- as.double(nU) * nX
  nMAX <- spatstat.options("maxmatrix")
  needsplit <- (nmat > nMAX)
  return(needsplit)
}

# function that should be called to evaluate interaction terms
# between quadrature points and data points

evalInteraction <- function(X, P, E = equalpairs(P, X), 
                            interaction, correction,
                            ...,
                            precomputed=NULL,
                            savecomputed=FALSE) {

  # evaluate the interaction potential
  # (does not assign/touch the variable names)

  verifyclass(interaction, "interact")

  # handle Poisson case
  if(is.poisson(interaction)) {
    out <- matrix(, nrow=npoints(P), ncol=0)
    attr(out, "IsOffset") <- logical(0)
    return(out)
  }
  
  # determine whether to use fast evaluation in C
  fastok    <- (spatstat.options("fasteval") %in% c("on", "test"))
  if(fastok) {
    cando   <- interaction$can.do.fast
    par     <- interaction$par
    dofast  <- !is.null(cando) && cando(X, correction, par)
  } else dofast <- FALSE

  # determine whether to split quadscheme into blocks
  if(dofast) {
    dosplit <- FALSE
  } else {
    # decide whether the quadrature scheme is too large to handle in one piece
    needsplit <- oversize.quad(nU=npoints(P), nX=npoints(X))

    # not implemented when savecomputed=TRUE
    dosplit   <- needsplit && !savecomputed
    if(needsplit && savecomputed)
      warning(paste("Oversize quadscheme cannot be split into blocks",
                    "because savecomputed=TRUE;",
                    "memory allocation error may occur"))
  }
  
  if(!dosplit) {
    # normal case
    V <- evalInterEngine(X=X, P=P, E=E,
                         interaction=interaction,
                         correction=correction,
                         ...,
                         precomputed=precomputed,
                         savecomputed=savecomputed)
  } else {
    # Too many quadrature points: split into blocks
    nX <- npoints(X)
    nP <- npoints(P)
    # Determine which evaluation points are data points
    Pdata <- E[,2]
    # hence which are dummy points
    Pall <- seq_len(nP)
    Pdummy <- if(length(Pdata) > 0) Pall[-Pdata] else Pall
    nD <- length(Pdummy)
    # size of full matrix
    nmat <- (nD + nX) * nX
    nMAX <- spatstat.options("maxmatrix")
    # Calculate number of dummy points in largest permissible X * (X+D) matrix 
    nperblock <- max(1, floor(nMAX/nX - nX))
    # determine number of such blocks 
    nblocks <- ceiling(nD/nperblock)
    nfull <- nblocks - 1
    # announce
    if(nblocks > 1) 
      message(paste("Large quadrature scheme",
                    "split into blocks to avoid memory size limits;",
                    nD, "dummy points",
                    "split into", nblocks, "blocks,",
                    "the first",
                    if(nfull > 1) paste(nfull, "blocks") else "block",
                    "containing",
                    nperblock, "dummy", ngettext(nperblock, "point", "points"),
                    "and the last block containing",
                    nD - nperblock * nfull, "dummy points"))
    #
    #
    Ei <- cbind(1:nX, 1:nX)
    #
    for(iblock in 1:nblocks) {
      first <- min(nD, (iblock - 1) * nperblock + 1)
      last  <- min(nD, iblock * nperblock)
      # extract dummy points  
      Di <- P[Pdummy[first:last]]
      Pi <- superimpose(X, Di, check=FALSE, W=X$window)
      # evaluate potential
      Vi <- evalInterEngine(X=X, P=Pi, E=Ei, 
                            interaction=interaction,
                            correction=correction,
                            ...,
                            savecomputed=FALSE)
      if(iblock == 1) {
        V <- Vi
      } else {
        # tack on the glm variables for the extra DUMMY points only
        V <- rbind(V, Vi[-(1:nX), , drop=FALSE])
      }
    }
  } 
  return(V)
}

# workhorse function that actually calls relevant code to evaluate interaction

evalInterEngine <- function(X, P, E, 
                            interaction, correction,
                            ...,
                            precomputed=NULL,
                            savecomputed=FALSE) {

  # fast evaluator (C code) may exist
  fasteval <- interaction$fasteval
  cando    <- interaction$can.do.fast
  par      <- interaction$par
  feopt    <- spatstat.options("fasteval")
  dofast   <- !is.null(fasteval) &&
              (is.null(cando) || cando(X, correction,par)) &&
              (feopt %in% c("on", "test"))
    
  V <- NULL
  if(dofast) {
    if(feopt == "test")
      message("Calling fasteval")
    V <- fasteval(X, P, E,
                  interaction$pot, interaction$par, correction, ...)
  }
  if(is.null(V)) {
    # use generic evaluator for family
    evaluate <- interaction$family$eval
    if("precomputed" %in% names(formals(evaluate))) {
      # Use precomputed data
      # version 1.9-3 onward (pairwise and pairsat families)
      V <- evaluate(X, P, E,
                    interaction$pot,
                    interaction$par,
                    correction, ...,
                    precomputed=precomputed,
                    savecomputed=savecomputed)
    } else {
      # Cannot use precomputed data
      # Object created by earlier version of ppm
      # or not pairwise/pairsat interaction
      V <- evaluate(X, P, E,
                    interaction$pot,
                    interaction$par,
                    correction, ...)
    }
  }

  return(V)
}


