#'
#'  rhohatlpp.R
#'
#'  rhohat.lpp and rhohat.lppm
#'
#'  Moved from rhohat.R to separate file rhohatlpp.R on 16 june 2020
#' 
#' Copyright (c) Adrian Baddeley 2015-2019
#' GNU Public Licence GPL >= 2.0

rhohat.lpp <- rhohat.lppm <- 
  function(object, covariate, ...,
           weights=NULL,
           method=c("ratio", "reweight", "transform"),
           horvitz=FALSE,
           smoother=c("kernel", "local", "decreasing", "increasing"),
           subset=NULL,
           nd=1000, eps=NULL, random=TRUE, 
           n=512, bw="nrd0", adjust=1, from=NULL, to=NULL, 
           bwref=bw, covname, confidence=0.95, positiveCI) {
  callstring <- short.deparse(sys.call())
  smoother <- match.arg(smoother)
  method <- match.arg(method)
  if(missing(positiveCI))
    positiveCI <- (smoother == "local")
  if(missing(covname)) 
    covname <- sensiblevarname(short.deparse(substitute(covariate)), "X")
  if(is.null(adjust))
    adjust <- 1
  # validate model
  if(is.lpp(object)) {
    X <- object
    model <- lppm(object, ~1, eps=eps, nd=nd, random=random)
    reference <- "Lebesgue"
    modelcall <- NULL
  } else if(inherits(object, "lppm")) {
    model <- object
    X <- model$X
    reference <- "model"
    modelcall <- model$call
  } else stop("object should be of class lpp or lppm")
  
  if("baseline" %in% names(list(...)))
    warning("Argument 'baseline' ignored: not available for ",
            if(is.lpp(object)) "rhohat.lpp" else "rhohat.lppm")

  if(is.character(covariate) && length(covariate) == 1) {
    covname <- covariate
    switch(covname,
           x={
             covariate <- function(x,y) { x }
           }, 
           y={
             covariate <- function(x,y) { y }
           },
           stop("Unrecognised covariate name")
         )
    covunits <- unitname(X)
  } else {
    covunits <- NULL
  }

  S <- as.psp(as.linnet(X))
  if(!is.null(subset)) S <- S[subset]
  totlen <- sum(lengths_psp(S))
  
  rhohatEngine(model, covariate, reference, totlen, ...,
               subset=subset,
               weights=weights,
               method=method,
               horvitz=horvitz,
               smoother=smoother,
               resolution=list(nd=nd, eps=eps, random=random),
               n=n, bw=bw, adjust=adjust, from=from, to=to,
               bwref=bwref, covname=covname, covunits=covunits,
               confidence=confidence, positiveCI=positiveCI,
               modelcall=modelcall, callstring=callstring)
}

