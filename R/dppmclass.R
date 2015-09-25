is.dppm <- function#Recognise Fitted Determinantal Point Process Models
### Check that an object inherits the class dppm
(x
 ### Any object.
 ){
    inherits(x, "dppm")
    ### A single logical value.

    ##keyword<< spatial
    ##keyword<< manip
    ##keyword<< models
}

plot.dppm <- function (x, ..., what = c("intensity", "statistic")){
    objectname <- short.deparse(substitute(x))
    if(missing(what) && is.stationary(x))
        what <- "statistic"
    plot.kppm(x, ..., objectname = objectname, what = what)
}

Kmodel.dppm <- function (model, ...){
    Kmodel(model$fitted, W=model$window)
}

pcfmodel.dppm <- function (model, ...){
    pcfmodel(model$fitted, W=model$window)
}

intensity.dppm <- function (X, ...){
    return(intensity(X$fitted))
}

reach.dppm <- function(x, ...){
    reach(x$fitted, ...)
}
