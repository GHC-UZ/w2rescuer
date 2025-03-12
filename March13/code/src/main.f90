program main
include 'param.h'

!---------opening of files for writing-------------
open(14,file='../output/xaxis.m')
open(11,file='../output/depth.m')
open(13,file='../output/bed.m')
open(15,file='../output/porosity.m')
open(33,file='../output/gauges.txt')
open(32,file='../output/surface.txt')
open(35,file='../output/mass.txt')
open(36,file='../output/energy.txt')



!-------read values from file incude.dat----------
open(3,file='../input/input.dat')
read(3,*) n
read(3,*) cfl
read(3,*) order
read(3,*) Tt
read(3,*) a
read(3,*) b
read(3,*) threshold
read(3,*) entropy
read(3,*) plt
read(3,*) limiter
read(3,*) Nm
read(3,*) theta
read(3,*) itime
read(3,*) amplitude
read(3,*) swl
close(3)

!------initialization----------------------
g=9.81d0
dx=dble((b-a)/(dble(n)-1.0d0))
x=a
t= 0.0d0
flag=0
pt=1
energyf(:)=0.0d0

call grid(dxi,problem,xc)
call initial(qp,n,dx,a,g,maxl,mass,etot,be,dxi,amplitude,swl,por)
call primitives(qp,pv,por,be)

t=0.0d0
ite=0
call popr(qp,be,dxi,g,t,ite,por) 

dt=(dx*cfl)/maxl
massin=mass
me = (massin-mass)/massin

call gauge(gp)

!---------Loop in time----------------------
if(itime.eq.1) then
        call euler(qp,pv,be,gp,ite,massin,por) 
else
        print*,'Only Eurer is implemented'
        stop
endif

!-----post processing----------------------- 
t=Tt
if(t.ne.0.0d0) call popr(qp,be,dxi,g,t,ite,por) 

stop
end
!*************************************************************************
double precision function pw(x)
double precision x

pw = dexp(dlog(x)/3.0d0)
return
end
!*************************************************************************
double precision function pow4(x)
double precision pw,x
pow4 = x*pw(x)
return
end
!*******************************************************************
subroutine primitives(qp,pv,por,bed)
implicit none
integer n,i
double precision dt,dx,a,b,L,eps_h,eps_u,Mass_dry,Vol_wet
double precision qp(4,-3:n+4),pv(4,-3:n+4),threshold
double precision por(-3:n+4),h_loc,hu_loc,bed(-3:n+4)
common /tr/ threshold
common /bd/ a,b
common /dn/ dt,dx,n
 
!-----------------version 1----------------------------
!do i=-3,n+4
!	pv(1,i)=qp(1,i)
!   	if(qp(1,i).lt.threshold) then 
!      		pv(2,i)=0.0d0
!      		pv(1,i)=0.0d0
!   	else
!      		pv(2,i)=qp(2,i)/qp(1,i)
!      		if(dabs(pv(2,i)).lt.1.d-6)pv(2,i)=0.0d0
!   	endif
!enddo
!
!-----------------version2------------------------
L = dabs(b-a)
eps_u = dx*dx/L
eps_h = 1.d-2*eps_u
threshold=eps_h
   
Mass_dry = 0.0d0
Vol_wet = 0.0d0
   
!do i=1,n
!	if (qp(1,i).lt.eps_h) then
!		Mass_dry = Mass_dry + qp(1,i) !total mass lost dry cells
!	else  
1		Vol_wet = Vol_wet + 1 ! total volume of wet areas
!	endif
!enddo

do i=-3,n+4
   	if (qp(1,i).lt.eps_h) then
      		qp(1,i) = 0
      		qp(2,i) = 0
      		pv(1,i) = 0
      		pv(2,i) = 0
   	else
!redistribute the lost mass
!      		qp(1,i) = qp(1,i) + (1.0/Vol_wet)*Mass_dry
!      		pv(1,i) = qp(1,i)

!**************** Ricchiuto approach *********************************
!      		if (qp(1,i).lt.eps_u) then
!         		qp(2,i) = 0
!         		pv(2,i) = 0
!      		else
!         		pv(2,i) = qp(2,i)/qp(1,i)
!      		endif
!************** Kourganov approach *******************************
                 h_loc = qp(1,i)/por(i) - bed(i)
                 hu_loc = qp(2,i)/por(i)
                 pv(1,i) = h_loc 
                 pv(2,i) = qp(2,i)/(qp(1,i)-qp(3,i)*qp(4,i)) !not corrected velocity/ problems for wet/dry fronts
      		 !pv(2,i) = h_loc*hc_loc/sqrt(max(h_loc**2,eps_u)**2)
      		 !qp(2,i) = qp(1,i)*pv(2,i)

   	endif
enddo

return  
end
!*************************************************************************
subroutine inflow(q,t)
implicit none
integer i,n
double precision q(2,-3:n+4),t,hta,dx,dt
double precision w,l,pi,am,k,dd,Tt,alpha,u1
double precision a,b,x,tmpam
common /dn/ dt,dx,n
common /bd/ a,b

pi=4.0d0*datan(1.0d0)
am=0.0215d0      
dd=0.36d0
l=6.141d0
Tt=3.3333d0
k=2.0d0*pi/l
alpha=-0.39d0
w=2.0d0*pi/Tt

x=a-4*dx
!do i=-3,1
tmpam=am
hta=tmpam*dsin(-w*t)
u1=hta*w/((1.0d0-(alpha+1.0d0/3.0d0)*(k*dd)**2)*k*dd)
q(1,1)=hta+dd
q(2,1)=q(1,1)*u1
!enddo

return
end
