#
# compareFit.R
#
# $Revision: 1.1 $  $Date: 2011/06/19 06:08:20 $

compareFit <- function(object, Fun, r=NULL, breaks=NULL,
                     ..., trend=~1, interaction=Poisson(),
                     rbord=NULL, modelnames=NULL,
                     same=NULL, different=NULL) {
  dotargs <- list(...)
  h <- hyperframe(obj=object, tren=trend, inte=interaction)
  N <- nrow(h)
  if(N == 0)
    stop("No objects specified")
  # determine rbord for summary statistics
  if(is.null(rbord) && !is.null(interaction))
    rbord <- max(with(h, reach(inte)))
  h$rbord <- rbord
  # try to get nice model names
  if(is.null(modelnames)) {
    if(inherits(trend, "formula") && is.interact(interaction) &&
       inherits(object, "listof") && all(nzchar(names(object))) &&
       length(names(object)) == nrow(h))
      modelnames <- names(object)
    else if(inherits(trend, "listof") && all(nzchar(names(trend))) &&
            length(names(trend)) == nrow(h))
      modelnames <- names(trend) 
    else if(inherits(interaction, "listof") &&
            all(nzchar(names(interaction))) &&
            length(names(interaction)) == nrow(h))
      modelnames <- names(interaction)
    else 
      modelnames <- row.names(h)
  }
  row.names(h) <- make.names(modelnames)
  # fix a common vector of r values
  if(is.null(r)) {
    # compute first function 
    fun1 <- with(h[1,,drop=FALSE],
                 do.call(Fun,
                         append(list(object=obj,
                                     trend=tren,
                                     interaction=inte,
                                     rbord=rbord,
                                     r=NULL, breaks=breaks),
                                dotargs)))
    # extract r values
    r <- with(fun1, .x)
  }
  # compute the subsequent functions
  if(N == 1)
    funs2toN <- NULL
  else 
    funs2toN <- with(h[-1, , drop=FALSE],
                     do.call(Fun,
                             append(list(object=obj,
                                         trend=tren,
                                         interaction=inte,
                                         rbord=rbord,
                                         r=r),
                                    dotargs)))
  if(N == 2)
    funs2toN <- list(funs2toN)
  # collect all functions in a list
  funs <- as.listof(append(list(fun1), funs2toN))
  names(funs) <- row.names(h)
  # collapse together
  out <- collapse.fv(funs, same=same, different=different)
  return(out)
}
