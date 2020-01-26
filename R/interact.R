#
#	interact.S
#
#
#	$Revision: 1.29 $	$Date: 2020/01/26 04:25:08 $
#
#	Class 'interact' representing the interpoint interaction
#               of a point process model
#              (e.g. Strauss process with a given threshold r)
#
#       Class 'isf' representing a generic interaction structure
#              (e.g. pairwise interactions)
#
#	These do NOT specify the "trend" part of the model,
#	only the "interaction" component.
#
#               The analogy is:
#
#                       glm()             ppm()
#
#                       model formula     trend formula
#
#                       family            interaction
#
#               That is, the 'systematic' trend part of a point process
#               model is specified by a 'trend' formula argument to ppm(),
#               and the interpoint interaction is specified as an 'interact'
#               object.
#
#       You only need to know about these classes if you want to
#       implement a new point process model.
#
#       THE DISTINCTION:
#       An object of class 'isf' describes an interaction structure
#       e.g. pairwise interaction, triple interaction,
#       pairwise-with-saturation, Dirichlet interaction.
#       Think of it as determining the "order" of interaction
#       but not the specific interaction potential function.
#
#       An object of class 'interact' completely defines the interpoint
#       interactions in a specific point process model, except for the
#       regular parameters of the interaction, which are to be estimated
#       by ppm() or otherwise. An 'interact' object specifies the values
#       of all the 'nuisance' or 'irregular' parameters. An example
#       is the Strauss process with a given, fixed threshold r
#       but with the parameters beta and gamma undetermined.
#
#       DETAILS:
#
#       An object of class 'isf' contains the following:
#
#	     $name               Name of the interaction structure         
#                                        e.g. "pairwise"
#
#	     $print		 How to 'print()' this object
#				 [A function; invoked by the 'print' method
#                                 'print.isf()']
#
#            $eval               A function which evaluates the canonical
#                                sufficient statistic for an interaction
#                                of this general class (e.g. any pairwise
#                                interaction.)
#
#       If lambda(u,X) denotes the conditional intensity at a point u
#       for the point pattern X, then we assume
#                  log lambda(u, X) = theta . S(u,X)
#       where theta is the vector of regular parameters,
#       and we call S(u,X) the sufficient statistic.
#
#       A typical calling sequence for the $eval function is
#
#            (f$eval)(X, U, E, potentials, potargs, correction)
#
#       where X is the data point pattern, U is the list of points u
#       at which the sufficient statistic S(u,X) is to be evaluated,
#       E is a logical matrix equivalent to (X[i] == U[j]),
#       $potentials defines the specific potential function(s) and
#       $potargs contains any nuisance/irregular parameters of these
#       potentials [the $potargs are passed to the $potentials without
#       needing to be understood by $eval.]
#       $correction is the name of the edge correction method.
#
#
#       An object of class 'interact' contains the following:
#
#
#            $name               Name of the specific potential
#                                        e.g. "Strauss"
#
#            $family              Object of class "isf" describing
#                                the interaction structure
#
#            $pot	         The interaction potential function(s)
#                                -- usually a function or list of functions.
#                                (passed as an argument to $family$eval)
#
#            $par                list of any nuisance/irregular parameters
#                                (passed as an argument to $family$eval)
#
#            $parnames           vector of long names/descriptions
#                                of the parameters in 'par'
#
#            $init()             initialisation action
#                                or NULL indicating none required
#
#            $update()           A function to modify $par
#                                [Invoked by 'update.interact()']
#                                or NULL indicating a default action
#
#	     $print		 How to 'print()' this object
#				 [Invoked by 'print' method 'print.interact()']
#                                or NULL indicating a default action
#
# --------------------------------------------------------------------------

print.isf <- function(x, ...) {
  if(is.null(x)) return(invisible(NULL))
  verifyclass(x, "isf")
  if(!is.null(x$print))
    (x$print)(x)
  invisible(NULL)
}

print.interact <- function(x, ..., family, brief=FALSE, banner=TRUE) {
  verifyclass(x, "interact")
  if(missing(family)) family <- waxlyrical('extras')
  #' Print name of model
  if(banner) {
    if(family && !brief && !is.null(xf <- x$family))
      print.isf(xf)
    splat(if(!brief) "Interaction:" else NULL, x$name, sep="")
  }
  # Now print the parameters
  if(!is.null(x$print)) {
     (x$print)(x)
  } else {
    # default
    # just print the parameter names and their values
    pwords <- x$parnames
    parval <- x$par
    pwords <- paste(toupper(substring(pwords, 1, 1)),
                    substring(pwords, 2), sep="")
    isnum <- sapply(parval, is.numeric)
    parval[isnum] <- lapply(parval[isnum], signif,
                            digits=getOption("digits"))
    splat(paste(paste0(pwords, ":\t", parval), collapse="\n"))
  }
  invisible(NULL)
}

is.interact <- function(x) { inherits(x, "interact") }

update.interact <- function(object, ...) {
  verifyclass(object, "interact")
  if(!is.null(object$update))
    (object$update)(object, ...)
  else {
    # Default
    # First update the version
    if(outdated.interact(object))
      object <- reincarnate.interact(object)
    # just match the arguments in "..."
    # with those in object$par and update them
    want <- list(...)
    if(length(want) > 0) {
      m <- match(names(want),names(object$par))
      nbg <- is.na(m)
      if(any(nbg)) {
        which <- paste((names(want))[nbg])
        warning(paste("Arguments not matched: ", which))
      }
      m <- m[!nbg]
      object$par[m] <- want
    }
    # call object's own initialisation routine
    if(!is.null(object$init))
      (object$init)(object)
    object
  }    
}

  
is.poisson.interact <- function(x) {
  verifyclass(x, "interact")
  is.null(x$family)
}

parameters.interact <- function(model, ...) {
  model$par
}

# Test whether interact object was made by an older version of spatstat

outdated.interact <- function(object) {
  ver <- object$version
  older <- is.null(ver) || (package_version(ver) < versionstring.spatstat())
  return(older)
}


# Test whether the functions in the interaction object
# expect the coefficient vector to contain ALL coefficients,
# or only the interaction coefficients.
# This change was introduced in 1.11-0, at the same time
# as interact objects were given version numbers.

newstyle.coeff.handling <- function(object) {
  stopifnot(inherits(object, "interact"))  
  ver <- object$version
  old <- is.null(ver) || (package_version(ver) < "1.11")
  return(!old)
}

# ######
#
# Re-create an interact object in the current version of spatstat
#
# 

reincarnate.interact <- function(object) {
  # re-creates an interact object in the current version of spatstat

  if(!is.null(object$update)) {
    newobject <- (object$update)(object)
    return(newobject)
  }
  
  par <- object$par
#  pot <- object$pot
  name <- object$name
  
  # get creator function
  creator <- object$creator
  if(is.null(creator)) {
    # old version: look up list
    creator <- .Spatstat.Old.InteractionList[[name]]
    if(is.null(creator))
      stop(paste("Don't know how to update", sQuote(name),
                 "to current version of spatstat"))
  }
  if(is.character(creator))
    creator <- get(creator)
  if(!is.function(creator) && !is.expression(creator))
    stop("Internal error: creator is not a function or expression")

  # call creator
  if(is.expression(creator)) 
    newobject <- eval(creator)
  else {
    # creator is a function
  
    # It's assumed that the creator function's arguments are
    # either identical to components of 'par' (the usual case)
    # or to one of the components of the object itself (Ord, Saturated)
    # or to printfun=object$print (Pairwise).
  
    argnames <- names(formals(creator))
    available <- append(par, object)
    available <- append(available, list(printfun=object$print))
    ok <- argnames %in% names(available)
    if(!all(ok))
      stop(paste("Internal error:",
                 ngettext(sum(!ok), "argument", "arguments"),
                 paste(sQuote(argnames[!ok]), collapse=", "),
                 "in creator function were not understood"))
    newobject <- do.call(creator, available[argnames])
  }
  
  if(!inherits(newobject, "interact"))
    stop("Internal error: creator did not return an object of class interact")

  return(newobject)
}


# This list is necessary to deal with older formats of 'interact' objects
# which did not include the creator name

.Spatstat.Old.InteractionList <-
  list("Diggle-Gratton process"    = "DiggleGratton",
       "Geyer saturation process"  = "Geyer",
       "Lennard-Jones potential"   = "LennardJones",
       "Multitype Strauss process" = "MultiStrauss",
       "Multitype Strauss Hardcore process" = "MultiStraussHard",
       "Ord process with threshold potential"="OrdThresh",
       "Piecewise constant pairwise interaction process"="PairPiece",
       "Poisson process"           = "Poisson",
       "Strauss process"           = "Strauss",
       "Strauss - hard core process" = "StraussHard",
       "Soft core process" = "Softcore",
       # weird ones:
       "Ord process with user-defined potential" = expression(Ord(object$pot)),
       "Saturated process with user-defined potential"
          =expression(Saturated(object$pot)),
       "user-defined pairwise interaction process"=
       expression(
           Pairwise(object$pot,
                    par=object$par,
                    parnames=object$parnames,
                    printfun=object$print))
           
     )
       
as.interact <- function(object) {
  UseMethod("as.interact")
}

as.interact.interact <- function(object) {
  verifyclass(object, "interact")
  return(object)
}

interactionfamilyname <- function(x) {
  if(inherits(x, "isf")) return(x$name)
  x <- as.interact(x)
  return(x$family$name)
}
                                      
#### internal code for streamlining initialisation of interactions
#
#    x should be a partially-completed 'interact' object
#

instantiate.interact <- function(x, par=NULL) {
  if(is.character(x$family)) x$family <- get(x$family)
  # set parameter values
  x$par    <- par
  # validate parameters
  if(!is.null(x$init)) x$init(x)
  x$version <- versionstring.spatstat()
  return(x)
}
