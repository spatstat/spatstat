#'
#'   bw.abram.R
#'
#'   Abramson bandwidths
#'
#'   $Revision: 1.3 $ $Date: 2019/03/05 06:41:39 $
#'

bw.abram <- function(X, h0, 
                     at=c("points", "pixels"), ...,
                     hp=NULL, pilot=NULL, trim=5){
  stopifnot(is.ppp(X))
  check.1.real(h0)
  check.1.real(trim)
  at <- match.arg(at)
  
  stopifnot(h0 > 0)
  stopifnot(trim > 0)
  
  pd <- pilot
  imwin <- as.im(Window(X), ...)
  pilot.data <- X
  
  if(!is.null(pd)){
    if(is.im(pd)){
      if(!compatible.im(imwin,pd))
        stop("'X' and 'pilot' have incompatible spatial domains", call.=FALSE)
      pilot[pd<=0] <- min(pd[pd>0]) # clip the worst small values away
      hp <- NULL
    } else if(is.ppp(pd)){
      if(!compatible.im(imwin,as.im(Window(pd), ...)))
        stop("'X' and 'pilot' have incompatible spatial domains", call.=FALSE)
      pilot.data <- pd
    } else {
      stop("if supplied, 'pilot' must be a pixel image or a point pattern",
           call.=FALSE)
    }
  }
  
  if(!is.im(pilot)){
    if(is.null(hp))
      stop(paste("pilot bandwidth 'hp' is required when 'pilot' is",
                 "unsupplied or given as point pattern"),
           call.=FALSE)
    pilot <- density(pilot.data,sigma=hp,positive=TRUE,...)
  }
  
  pilot <- pilot/integral(pilot) # scale to probability density
  pilotvalues <- safelookup(pilot, pilot.data, warn=FALSE)
  ## geometric mean re-scaler (Silverman, 1986; ch 5).  
  gamma <- exp(mean(log(pilotvalues[pilotvalues > 0])))^(-0.5)

  switch(at,
         points = {
           pilot.X <- safelookup(pilot,X,warn=FALSE)
           bw <- h0 * pmin((pilot.X^(-0.5))/gamma,trim)
         },
         pixels = {
           bw <- eval.im(h0 * pmin((pilot^(-0.5))/gamma, trim))
         })
  
  return(bw)
}
