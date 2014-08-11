C Output from Public domain Ratfor, version 1.0
      subroutine inxyp(x,y,xp,yp,npts,nedges,score,onbndry)
      implicit double precision(a-h,o-z)
      dimension x(npts), y(npts), xp(nedges), yp(nedges), score(npts)
      logical first, onbndry(npts)
      zero = 0.0d0
      half = 0.5d0
      one = 1.0d0
      do23000 i = 1,nedges 
      
      x0 = xp(i)
      y0 = yp(i)
      if(i .eq. nedges)then
      
      x1 = xp(1)
      y1 = yp(1)
      else
      
      x1 = xp(i+1)
      y1 = yp(i+1)
      endif
      
      dx = x1 - x0
      dy = y1 - y0
      do23004 j = 1,npts 
      
      xcrit = (x(j) - x0)*(x(j) - x1)
      if(xcrit .le. zero)then
      
      if(xcrit .eq. zero)then
      
      contrib = half
      else
      
      contrib = one
      endif
      
      ycrit = y(j)*dx - x(j)*dy + x0*dy - y0*dx
      if(dx .lt. 0)then
      
      if(ycrit .ge. zero)then
      
      score(j) = score(j) + contrib
      endif
      
      onbndry(j) = onbndry(j) .or. (ycrit .eq. zero)
      else
      if(dx .gt. zero)then
      
      if(ycrit .lt. zero)then
      
      score(j) = score(j) - contrib
      endif
      
      onbndry(j) = onbndry(j) .or. (ycrit .eq. zero)
      else
      
      if(x(j) .eq. x0)then
      
      ycrit = (y(j) - y0)*(y(j) - y1)
      endif
      
      onbndry(j) = onbndry(j) .or. (ycrit .le. zero)
      endif
      endif
      
      endif
      
23004 continue
23005 continue
      
23000 continue
23001 continue
      return
      end
