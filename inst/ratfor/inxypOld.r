subroutine inxyp(x,y,xp,yp,npts,nedges,score,onbndry)
implicit double precision(a-h,o-z)
dimension x(npts), y(npts), xp(nedges), yp(nedges), score(npts)
logical first, onbndry(npts)
zero = 0.0d0
half = 0.5d0
one  = 1.0d0
do i = 1,nedges {
  x0 = xp(i)
  y0 = yp(i)
  if(i == nedges) {
    x1 = xp(1)
    y1 = yp(1)
  } else {
    x1 = xp(i+1)
    y1 = yp(i+1)
  }
  dx = x1 - x0
  dy = y1 - y0
  do j = 1,npts {
    xcrit = (x(j) - x0)*(x(j) - x1)
    if(xcrit <= zero) {
      if(xcrit == zero) {
        contrib = half
      } else {
        contrib = one
      }
      ycrit = y(j)*dx - x(j)*dy + x0*dy - y0*dx
      if(dx < 0) {
        if(ycrit >= zero) {
          score(j) = score(j) + contrib
        }
        onbndry(j) = onbndry(j) | (ycrit == zero)
      } else if(dx > zero) {
        if(ycrit < zero) {
          score(j) = score(j) - contrib
        }
        onbndry(j) = onbndry(j) | (ycrit == zero)
      } else {
        if(x(j) == x0) {
          ycrit = (y(j) - y0)*(y(j) - y1)
        }
        onbndry(j) = onbndry(j) | (ycrit <= zero)
      }
    }
  }
}
return
end
