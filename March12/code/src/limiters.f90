double precision function fhi(fh,cl,cr,x,limiter,bb)
implicit none
integer n,limiter,tmplim,bb(n),cl,cr
double precision minmod
double precision x,y,dt,dx,fh(0:n)
common /dn/ dt,dx,n


if(limiter.eq.1) then
	fhi = minmod(1.0d0,x)
elseif(limiter.eq.2) then
  	fhi = dmax1(0.0d0,dmin1(1.0d0,2.0d0*x),dmin1(2.0d0,x))
elseif(limiter.eq.3) then
  	fhi = dmax1(0.0d0,dmin1((1.0d0+x)/2.0d0,2.0d0,2.0d0*x))
elseif(limiter.eq.4) then
  	fhi = (x + dabs(x))/(1.0d0+dabs(x))
elseif(limiter.eq.5) then
  	!fhi = (2.0d0*x)/(1.0d0+x**2) !van albada
  	fhi = (x**2+x)/(1.0d0+x**2) !van albada
elseif(limiter.eq.6) then
  	if(x.le.0) then
     		fhi=0.0d0
  	elseif( x.gt.0.0d0 .and. x.le.0.5d0) then 
     		fhi=2.0d0*x
  	elseif( x.gt.0.5d0 .and. x.le.1.0d0) then
     		fhi=1.0d0
  	elseif( x.gt.1.0d0) then
     		fhi=dmin1(x,2.0d0,2.0d0/(1.0d0+x))
  	endif
elseif(limiter.eq.7) then
 	 fhi=1.0d0
endif

return
end
!------------------------------------------------------  
double precision function minmod3(a,b,c)
implicit none
double precision a,b,c,s
s=dsign(1.0d0,a)
minmod3=s*dmax1(0.0d0,dmin1(dabs(a),2.0d0*s*b,2.0d0*s*c))
return
end
!------------------------------------------------------  
double precision function minmod(a,b)
implicit none
double precision a,b

if(dabs(a).lt.dabs(b) .and. (a*b).gt.0.0d0) then
   	minmod = a
elseif(dabs(b).lt.dabs(a) .and. (a*b).gt.0.0d0) then
   	minmod = b
elseif((a*b).le.0.0d0) then
   	minmod = 0.0d0
endif

return
end
