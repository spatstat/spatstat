#
#
#   rmhResolveTypes.R
#
#   $Revision: 1.10 $   $Date: 2019/02/20 03:34:50 $
#
#
rmhResolveTypes <- function(model, start, control) {

# Decide whether a multitype point process is to be simulated.
# If so, determine the vector of types.

  verifyclass(model, "rmhmodel")
  verifyclass(start, "rmhstart")
  verifyclass(control, "rmhcontrol")

# Different ways of specifying types directly

  types.model <- model$types
  types.start <- if(start$given=="x" && is.marked(x.start <- start$x.start))
                     levels(marks(x.start, dfok=FALSE)) else NULL
  
# Check for inconsistencies  
  if(!is.null(types.model) && !is.null(types.start))
    if(!isTRUE(all.equal(types.model, types.start)))
      stop("marks in start$x.start do not match model$types")
  
  types.given <- if(!is.null(types.model)) types.model else types.start
  types.given.source <-
    if(!is.null(types.model)) "model$types" else "marks of x.start"
  
# Different ways of implying the number of types
  
  ntypes.beta <- length(model$par[["beta"]])
  ntypes.ptypes <- length(control$ptypes)
  ntypes.nstart <- if(start$given == "n") length(start$n.start) else 0
  mot <- model$trend
  ntypes.trend <-  if(is.null(mot)) 0 else
                   if(is.im(mot)) 1 else
                   if(is.list(mot) &&
                      all(unlist(lapply(mot, is.im))))
                     length(mot) else 0
  
# Check for inconsistencies in implied number of types (only for numbers > 1)

  nty <- c(ntypes.beta, ntypes.ptypes, ntypes.nstart, ntypes.trend)
  nam <- c("model$par$beta", "control$ptypes", "start$n.start", "model$trend")
  implied <- (nty > 1)
  if(!any(implied))
    ntypes.implied <- 1
  else {
    if(length(unique(nty[implied])) > 1)
      stop(paste("Mismatch in numbers of types implied by",
                 commasep(sQuote(nam[implied]))))
    ntypes.implied <- unique(nty[implied])
    ntypes.implied.source <- (nam[implied])[1]
  } 

# Check consistency between types.given and ntypes.implied 

  if(!is.null(types.given) && ntypes.implied > 1)
    if(length(types.given) != ntypes.implied)
      stop(paste("Mismatch between number of types in",
                 types.given.source,
                 "and length of",
                 ntypes.implied.source))

# Finally determine the types
  
  if(model$multitype.interact) {
    # There MUST be a types vector
    types <- if(!is.null(types.given)) types.given
             else if(ntypes.implied > 1) 1:ntypes.implied
             else stop("Cannot determine types for multitype process")
  } else {
    types <- if(!is.null(types.given)) types.given
             else if(ntypes.implied > 1) 1:ntypes.implied
             else 1
  }

  ntypes <- length(types)
  
# If we are conditioning on the number of points of each type,
# make sure the starting state is appropriate

  if(control$fixing == "n.each.type") {
    if(start$given == "n" && ntypes.nstart != ntypes)
      stop("Length of start$n.start not equal to number of types.\n")
    else if(start$given == "x" && length(types.given) != ntypes) 
      stop("Marks of start$x.start do not match number of types.\n")
  }
  
  return(types)
}

  
