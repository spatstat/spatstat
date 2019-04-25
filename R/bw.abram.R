#'
#'   bw.abram.R
#'
#'   Abramson bandwidths
#'
#'   $Revision: 1.8 $ $Date: 2019/04/25 05:34:49 $
#'

bw.abram <- function(X, h0,
                     ..., 
                     at=c("points", "pixels"), 
                     hp=h0, pilot=NULL, trim=5,
                     smoother=density.ppp){
  stopifnot(is.ppp(X))
  at <- match.arg(at)

  if(missing(h0) || is.null(h0)) {
    h0 <- bw.ppl(X)
  } else {
    check.1.real(h0)
    stopifnot(h0 > 0)
  }

  check.1.real(trim)
  stopifnot(trim > 0)
  
  pilot.data <- X
  imwin <- as.im(Window(X), ...)
  
  if(is.im(pilot)){
    if(!compatible.im(imwin,pilot))
      stop("'X' and 'pilot' have incompatible spatial domains", call.=FALSE)
    #' clip the worst small values away
    pilot[pilot<=0] <- min(pilot[pilot>0]) 
  } else if(is.ppp(pilot)){
    if(!compatible.im(imwin,as.im(Window(pilot), ...)))
      stop("'X' and 'pilot' have incompatible spatial domains", call.=FALSE)
    pilot.data <- pilot
  } else if(!is.null(pilot))
    stop("if supplied, 'pilot' must be a pixel image or a point pattern",
         call.=FALSE)
  
  if(!is.im(pilot)) {
    if(is.character(smoother)) {
      smoother <- get(smoother, mode="function")
    } else stopifnot(is.function(smoother))
    pilot <- smoother(pilot.data,sigma=hp,positive=TRUE,...)
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
