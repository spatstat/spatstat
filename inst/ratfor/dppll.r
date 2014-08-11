subroutine dppll(x,y,l1,l2,l3,l4,np,nl,eps,mint,rslt,xmin,jmin)
implicit double precision(a-h,o-z)
dimension x(np), y(np), rslt(np,nl), xmin(np), jmin(np)
double precision l1(nl), l2(nl), l3(nl), l4(nl)
one = 1.d0
zero = 0.d0
do j = 1,nl {
	dx = l3(j) - l1(j)
	dy = l4(j) - l2(j)
	alen = sqrt(dx**2 + dy**2)
	if(alen .gt. eps) {
		co = dx/alen
		si = dy/alen
	} else {
          co = 0.5
          si = 0.5
        }
	do  i = 1, np {
		xpx1 = x(i) - l1(j)
		ypy1 = y(i) - l2(j)
		xpx2 = x(i) - l3(j)
		ypy2 = y(i) - l4(j)
		d1 = xpx1**2 + ypy1**2
		d2 = xpx2**2 + ypy2**2
		dd = min(d1,d2)
		if(alen .gt. eps) {
			xpr = xpx1*co + ypy1*si
			if(xpr .lt. zero .or. xpr .gt. alen) {
				d3 = -one
			}
			else {
				ypr = - xpx1*si + ypy1*co
				d3 = ypr**2
			}
		}
		else {
				d3 = -one
		}
		if(d3 .ge. zero) {
			dd = min(dd,d3)
		}
		sd =sqrt(dd)
		rslt(i,j) = sd
		if(mint.gt.0) {
			if(sd .lt. xmin(i)) {
				xmin(i) = sd
				if(mint.gt.1) {
					jmin(i) = j
				}
			}
		}
	}
}
return
end
