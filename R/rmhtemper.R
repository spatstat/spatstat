#'
#'      rmhtemper.R
#'
#'   $Revision: 1.2 $  $Date: 2015/08/29 04:36:38 $
#'

reheat <- local({

  expon <- function(x, alpha) {
    if(is.null(x)) return(NULL)
    if(is.numeric(x)) return(x^alpha)
    if(is.im(x)) return(x^alpha)
    if(is.function(x)) {
      f <- x
      g <- function(...) { f(...)^alpha }
      if(!inherits(f, "funxy")) return(g)
      return(funxy(g, W=as.owin(f)))
    }
    if(is.list(x)) return(lapply(x, expon))
    stop("Unrecognised format for x in x^alpha", call.=FALSE)
  }
    
  reheat <- function(model, invtemp) {
    model <- rmhmodel(model)
    cif   <- model$cif
    par   <- model$par
    w     <- model$w
    trend <- model$trend
    types <- model$types

    newtrend <- expon(trend, invtemp)

    rules <- lapply(cif, spatstatRmhInfo)
    temperfuns <- lapply(rules, getElement, name="temper")
    if(any(bad <- sapply(temperfuns, is.null)))
      stop(paste("reheating the", commasep(sQuote(cif[bad])),
                 ngettext(sum(bad), "cif", "cifs"),
                 "is not supported"))

    Ncif <- length(cif)
    if(Ncif == 1) {
      newpar <- temperfuns[[1]](par, invtemp)
    } else {
      newpar <- par
      for(i in 1:Ncif) 
        newpar[[i]] <- temperfuns[[i]](par[[i]], invtemp)
    }
    newmodel <- rmhmodel(cif=cif,
                         par=newpar, trend=newtrend,
                         w=w, types=types)
    return(newmodel)
  } 

  reheat
  
})


rtemper <- function(model, invtemp, nrep, ..., start=NULL, verbose=FALSE){
  df <- data.frame(invtemp, nrep)
  ndf <- nrow(df)
  X <- NULL
  for(i in 1:ndf) {
    if(verbose)
      cat(paste("Running", nrep[i], "steps",
                "at inverse temperature", invtemp[i], "... "))
    model.i <- reheat(model, invtemp[i])
    X <- rmh(model.i, nrep=nrep[i], ...,
             start=start,
             overrideXstart = X,
             overrideclip   = (i != ndf),
             verbose=FALSE)
    if(verbose) cat(" Done.\n")
  }
  return(X)
}
