##
## convexify.R
##
## $Revision: 1.2 $ $Date: 2020/03/16 10:28:51 $

convexify <- function(W, eps) {
  if(!is.polygonal(W)) {
    if(missing(eps)) eps <- diameter(Frame(W))/20
    W <- simplify.owin(W, eps)
  }
  e <- edges(W)
  len <- lengths_psp(e)
  ang <- angles.psp(e, directed=TRUE)
  df <- data.frame(ang=ang, len=len)
  df <- df[order(df$ang), ]
  df <- within(df, { dx <- len * cos(ang); dy <- len * sin(ang)})
  owin(poly=with(df, list(x=cumsum(c(0,dx)), y=cumsum(c(0,dy)))))
}

    
