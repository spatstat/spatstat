C Output from Public domain Ratfor, version 1.0
      subroutine dppll(x,y,l1,l2,l3,l4,np,nl,eps,mint,rslt,xmin,jmin)
      implicit double precision(a-h,o-z)
      dimension x(1), y(1), rslt(np,1), xmin(1), jmin(1)
      double precision l1(1), l2(1), l3(1), l4(1)
      one = 1.d0
      zero = 0.d0
      do23000 j = 1,nl 
      dx = l3(j) - l1(j)
      dy = l4(j) - l2(j)
      alen = sqrt(dx**2 + dy**2)
      if(alen .gt. eps)then
      co = dx/alen
      si = dy/alen
      else
      co = 0.5
      si = 0.5
      endif
      do23004 i = 1, np 
      xpx1 = x(i) - l1(j)
      ypy1 = y(i) - l2(j)
      xpx2 = x(i) - l3(j)
      ypy2 = y(i) - l4(j)
      d1 = xpx1**2 + ypy1**2
      d2 = xpx2**2 + ypy2**2
      dd = min(d1,d2)
      if(alen .gt. eps)then
      xpr = xpx1*co + ypy1*si
      if(xpr .lt. zero .or. xpr .gt. alen)then
      d3 = -one
      else
      ypr = - xpx1*si + ypy1*co
      d3 = ypr**2
      endif
      else
      d3 = -one
      endif
      if(d3 .ge. zero)then
      dd = min(dd,d3)
      endif
      sd =sqrt(dd)
      rslt(i,j) = sd
      if(mint.gt.0)then
      if(sd .lt. xmin(i))then
      xmin(i) = sd
      if(mint.gt.1)then
      jmin(i) = j
      endif
      endif
      endif
23004 continue
23005 continue
23000 continue
23001 continue
      return
      end
