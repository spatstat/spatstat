##
## convexify.R
##
## $Revision: 1.1 $ $Date: 2015/10/23 12:34:17 $

convexify <- function(W, eps) {
  if(!is.polygonal(W)) {
    if(missing(eps)) eps <- diameter(Frame(W))/20
    W <- simplify.owin(W, eps)
  }
  e <- edges(W)
  len <- lengths.psp(e)
  ang <- angles.psp(e, directed=TRUE)
  df <- data.frame(ang=ang, len=len)
  df <- df[order(df$ang), ]
  df <- within(df, { dx <- len * cos(ang); dy <- len * sin(ang)})
  owin(poly=with(df, list(x=cumsum(c(0,dx)), y=cumsum(c(0,dy)))))
}

    
