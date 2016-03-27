#'
#' rppm.R
#' 
#'  Recursive Partitioning for Point Process Models
#'
#'  $Revision: 1.7 $  $Date: 2016/03/27 10:18:07 $

rppm <- function(..., rpargs=list()) {
  fit <- ppm(..., forcefit=TRUE)
  if(!is.poisson(fit))
    warning("Interpoint interaction will be ignored", call.=FALSE)
  df <- getglmdata(fit)
  gf <- getglmfit(fit)
  rp <- do.call(rpart,
                resolve.defaults(list(formula=formula(gf), data=df),
                                 rpargs,
                                 list(method="poisson")))
  result <- list(fit=fit, rp=rp)
  class(result) <- c("rppm", class(result))
  return(result)
}

print.rppm <- function(x, ...) {
  splat("Point process model with recursive partitioning")
  splat("Data:", sQuote(x$fit$Qname))
  splat("Covariates:", commasep(sQuote(variablesinformula(formula(x$fit)))))
  splat("Regression tree:")
  print(x$rp)
  invisible(NULL)
}

plot.rppm <- function(x, ..., what=c("tree", "spatial")) {
  xname <- short.deparse(substitute(x))
  what <- match.arg(what)
  switch(what,
         tree = {
           out <- plot(x$rp, ...)
           text(x$rp, ...)
         },
         spatial = {
           p <- predict(x)
           out <- do.call("plot",
                          resolve.defaults(list(x=p),
                                           list(...),
                                           list(main=xname)))
         })
  return(invisible(out))
}

#' prune method

prune.rppm <- function(tree, ...) {
  tree$rp <- rpart::prune(tree$rp, ...)
  return(tree)
}

#' predict method

predict.rppm <- function(object, ...) {
  model <- object$fit
  tree  <- object$rp
  #' assemble covariates for prediction, using rules of predict.ppm
  co <- predict(model, ..., type="covariates", check=FALSE, repair=FALSE)
  newdata <- co$newdata
  masque  <- co$mask
  #' perform prediction using the tree
  pred <- predict(tree, newdata=co$newdata)
  #' pack up appropriately
  if(is.null(masque))
    return(pred)
  imago <- as.im(masque, value=1.0)
  if(!is.marked(model)) {
    out <- imago
    out[] <- pred
  } else {
    lev <- levels(marks(data.ppm(model)))
    nlev <- length(lev)
    out <- rep(list(imago), nlev)
    names(out) <- lev
    splitpred <- split(pred, newdata$marks)
    for(i in seq_len(nlev))
      out[[i]][] <- splitpred[[i]]
    out <- as.solist(out)
  }
  return(out)
}
    
fitted.rppm <- function(object, ...) {
  predict(object, locations=data.ppm(object$fit))
}

