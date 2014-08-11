#  Variance-covariance matrix for mppm objects
#
# $Revision: 1.9 $ $Date: 2013/11/09 16:55:46 $
#
#
vcov.mppm <- function(object, ..., what="vcov", err="fatal") {

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
    
  whinge <- NULL
  if(object$Fit$fitter != "glm")
    whinge <- "vcov.mppm only implemented for glm fits"
  else if(!is.poisson.mppm(object))
    whinge <- "vcov.mppm only implemented for Poisson processes"

  if(!is.null(whinge))
    return(errhandler(whinge, err))
  
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
