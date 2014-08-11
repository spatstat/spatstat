#
# fii.R
#
# Class of fitted interpoint interactions
#
#
fii <- function(interaction=NULL, coefs=numeric(0),
                Vnames=character(0), IsOffset=NULL) {
  if(is.null(interaction)) 
    interaction <- Poisson()
  stopifnot(is.interact(interaction))
  if(is.poisson.interact(interaction)) {
    if(length(Vnames) > 0)
      stop("Coefficients inappropriate for Poisson process")
  }
  if(is.null(IsOffset))
    IsOffset <- rep.int(FALSE, length(Vnames))
  else {
    stopifnot(is.logical(IsOffset))
    stopifnot(length(IsOffset) == length(Vnames))
  } 
  out <- list(interaction=interaction,
              coefs=coefs,
              Vnames=Vnames,
              IsOffset=IsOffset)
  class(out) <- c("fii", class(out))
  return(out)
}

summary.fii <- function(object, ...) {
  y <- unclass(object)
  INTERACT <- object$interaction
  coefs    <- object$coefs
  Vnames   <- object$Vnames
  IsOffset <- object$IsOffset
  y$poisson <- is.poisson.interact(INTERACT)
  thumbnail <- NULL
  if(y$poisson) {
    thumbnail <- "Poisson()"
  } else {
    if(!is.null(INTERACT$interpret)) {
      # invoke auto-interpretation feature
      sensible <-  
        if(newstyle.coeff.handling(INTERACT))
          (INTERACT$interpret)(coefs[Vnames[!IsOffset]], INTERACT)
        else 
          (INTERACT$interpret)(coefs, INTERACT)
      if(!is.null(sensible)) {
        header <- paste("Fitted", sensible$inames)
        printable <- sensible$printable
        # Try to make a thumbnail description
        param <- sensible$param
        ipar <- INTERACT$par
        if(all(unlist(lapply(param, length)) == 1) &&
           all(unlist(lapply(ipar, length)) == 1)) {
          allargs <- append(ipar, param)
          allargs <- lapply(allargs, signif, digits=4)
          thumbnail <- fakecallstring(INTERACT$creator, allargs)
        } 
      } else {
        # no fitted interaction parameters (e.g. Hard Core)
        header <- NULL
        printable <- NULL
        thumbnail <- paste0(INTERACT$creator, "()")
      }
    } else {
      # fallback
      sensible <- NULL
      VN <- Vnames[!IsOffset]
      if(length(VN) > 0) {
        header <- "Fitted interaction terms"
        icoef <- coefs[VN]
        printable <-  exp(unlist(icoef))
        ricoef <- lapply(icoef, signif, digits=4)
        thumbnail <- fakecallstring(INTERACT$creator, ricoef)
      } else {
        header <- NULL
        printable <- NULL
        thumbnail <- paste0(INTERACT$creator, "()")
      }
    }
    y <- append(y, list(sensible=sensible,
                        header=header,
                        printable=printable,
                        thumbnail=thumbnail))
  }
  class(y) <- c("summary.fii", class(y))
  return(y)
}

print.fii <- function(x, ...) {
  tiny <- resolve.1.default("tiny", list(...), list(tiny=FALSE))
  print(summary(x), brief=TRUE, tiny=tiny)
  return(invisible(NULL))
}

print.summary.fii <- function(x, ...) {
  secret <- resolve.defaults(list(...),
                             list(prefix="Interaction: ",
                                  family=TRUE,
                                  brief=FALSE,
                                  tiny=FALSE))
  if(secret$tiny) {
    # use thumbnail if available
    thumbnail <- x$thumbnail
    if(!is.null(thumbnail)) {
      cat(paste(thumbnail, "\n"))
      return(invisible(NULL))
    }
  }
  brief <- secret$brief
  if(!brief)
    cat(secret$prefix)
  if(x$poisson)
    cat("Poisson process\n")
  else {
    print(x$interaction, family=secret$family, brief=TRUE)
    if(!is.null(x$printable)) {
      nvalues <- length(x$printable)
      nheader <- length(x$header)
      if(nvalues == 1) {
        cat(paste(x$header, ":\t", x$printable, "\n", sep=""))
      } else if(nvalues == nheader) {
        for(i in 1:nheader) {
          cat(x$header[i])
          xpi <- x$printable[[i]]
          if(!is.list(xpi) && length(xpi) == 1) {
            cat(":\t", xpi, "\n")
          } else {
            cat(":\n")
            print(xpi)
          }
        } 
      } else {
        cat(paste(x$header, ":\n", sep=""))
        print(x$printable)
      } 
    }
  }
  if(!brief) {
    co <- x$coefs[x$Vnames[!x$IsOffset]]
    if(length(co) > 0) {
      cat("\nRelevant coefficients:\n")
      print(co)
    }
  }
  return(invisible(NULL))
}

coef.summary.fii <- function(object, ...) {
  object$printable
}

reach.fii <- function(x, ..., epsilon=0) {
  inte <- x$interaction
  coeffs <- x$coefs
  Vnames <- x$Vnames

  if(is.poisson.interact(inte))
    return(0)

  # get 'irange' function from interaction object
  irange <- inte$irange

  if(is.null(irange))
    return(Inf)

  # apply 'irange' function using fitted coefficients
  if(newstyle.coeff.handling(inte))
    ir <- irange(inte, coeffs[Vnames], epsilon=epsilon)
  else 
    ir <- irange(inte, coeffs, epsilon=epsilon)
  
  if(is.na(ir))
    ir <- Inf

  return(ir)
}

plot.fii <- function(x, ...) {
  if(is.poisson.interact(x$interaction)) {
    message("Poisson interaction; nothing plotted")
    return(invisible(NULL))
  }
  plfun <- x$interaction$family$plot
  if(is.null(plfun)) 
    stop("Plotting not implemented for this type of interaction")
  plfun(x, ...)
}


fitin <- function(object) {
  UseMethod("fitin")
}

fitin.ppm <- function(object) {
  f <- object$fitin
  if(!is.null(f))
    return(f)
  # For compatibility with older versions
  inte <- object$interaction
  if(is.null(inte)) 
    f <- fii() # Poisson
  else {
    coefs <- coef(object)
    Vnames <- object$internal$Vnames
    IsOffset <- object$internal$IsOffset
    # Internal names of regressor variables 
    f <- fii(inte, coefs, Vnames, IsOffset)
  }
  return(f)
}

as.interact.fii <- function(object) {
  verifyclass(object, "fii")
  return(object$interaction)
}

coef.fii <- function(object, ...) {
  verifyclass(object, "fii")
  return(object$coefs)
}
