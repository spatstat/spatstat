#
#    rescue.rectangle.R
# 
#    $Revision: 1.7 $   $Date: 2018/09/28 05:05:49 $
#
rescue.rectangle <- function(W) {
  verifyclass(W, "owin")

  if(W$type == "mask" && all(W$m))
     return(owin(W$xrange, W$yrange, unitname=unitname(W)))

  if(W$type == "polygonal" && length(W$bdry) == 1) {
    x <- W$bdry[[1]]$x
    y <- W$bdry[[1]]$y
    if(length(x) == 4 && length(y) == 4) {
      # could be a rectangle
      ux <- veryunique(x)
      uy <- veryunique(y)
      if(length(ux) == 2 && length(uy) == 2)
        return(owin(ux,uy, unitname=unitname(W)))
    }
  }
  return(W)
}

veryunique <- function(z) {
  uz <- sortunique(z)
  epsilon <- 2 * .Machine$double.eps * diff(range(uz))
  close <- (diff(uz) <= epsilon)
  uz <- uz[c(TRUE, !close)]
  return(uz)
}

