subroutine rk3(qp,pv,be,gp,ite,massin)
include 'param.h'
integer TSTEPS
double precision ARK(3,3),tmp_phi(n),PLTT(13)

!---------SYNOLAKIS---------------------
  PLTT(1)=dsqrt(1.0d0/9.81d0)*10.0d0
  PLTT(2)=dsqrt(1.0d0/9.81d0)*15.0d0
  PLTT(3)=dsqrt(1.0d0/9.81d0)*20.0d0
  PLTT(4)=dsqrt(1.0d0/9.81d0)*25.0d0
  PLTT(5)=dsqrt(1.0d0/9.81d0)*30.0d0
  PLTT(6)=dsqrt(1.0d0/9.81d0)*45.0d0
  PLTT(7)=dsqrt(1.0d0/9.81d0)*55.0d0
  PLTT(8)=dsqrt(1.0d0/9.81d0)*60.0d0
  PLTT(9)=dsqrt(1.0d0/9.81d0)*70.0d0
  PLTT(10)=dsqrt(1.0d0/9.81d0)*80.0d0

TSTEPS=3
ARK(1,1) = 1.0d0
ARK(1,2) = 0.0d0
ARK(1,3) = 1.0d0
ARK(2,1) = 3.0d0/4.0d0
ARK(2,2) = 1.0d0/4.0d0
ARK(2,3) = 1.0d0/4.0d0
ARK(3,1) = 1.0d0/3.0d0
ARK(3,2) = 2.0d0/3.0d0
ARK(3,3) = 2.0d0/3.0d0
 
if(Tt.eq.0.0d0) goto 110
g=9.81d0
iter=1
itert=0.0d0
ITE=0
it=0
pt=1
Q(1,:)=QP(1,:)
Q(2,:)=QP(2,:)

100  t = t + dt     

flagtmp=1
mass=0.0d0
etot=0.0d0

print*,'XRONOS: ' ,t , ' dt: ',dt    
ITE=ITE+1
!*************************************************************************

x=a
flag=0
maxl=0.0d0
helpmax=0.0d0
maxv=0.0d0

DO I=-3,N+4
	QP(1,I)=Q(1,I)
   	QP(2,I)=Q(2,I)
   	QSTAR(1,I)=Q(1,I)
   	QSTAR(2,I)=Q(2,I)
ENDDO

DO IRK=1,TSTEPS

	DO I=-3,N+4
   		QSTAR(1,I)=Q(1,I)
   		QSTAR(2,I)=Q(2,I)
	ENDDO

!---------everything is computed using values from the previous time step
 	call primitives(qstar,pv)
 	call flx_roe(qstar,pv,g,fl,entropy,order,limiter,sh,be,bb,fh,nutot,energyf) !calculate the fluxes ADVECTIVE PART
 	x=a
 	do i=1,n
    		F(1,i) = -(1.0d0/DX)*(fl(1,i)-fl(1,i-1)-sh(1,i)) 
    		F(2,i) = -(1.0d0/DX)*(fl(2,i)-fl(2,i-1)-sh(2,i))
 		x=x+dx
 	enddo

	DO I=1,N
   		Q(1,I) = ARK(IRK,1)*QP(1,I)+ARK(IRK,2)*QSTAR(1,I)+ARK(IRK,3)*DT*F(1,I)
   		Q(2,I) =  ARK(IRK,1)*QP(2,I)+ARK(IRK,2)*QSTAR(2,I)+ARK(IRK,3)*DT*F(2,I)
 	ENDDO

	call friction(q,pv,g,fr)

	call boundarya(q,n,t)
 	call boundaryb(q,n)



	DO I=-3,N+4
     		QSTAR(1,I)=Q(1,I)
     		QSTAR(2,I)=Q(2,I)
 	ENDDO 
    
ENDDO

!      call runup(q,be,n,t)
!************************************************************************
call primitives(q,pv)
DO I=1,N
	if(dsqrt(pv(1,i))+dsqrt(pv(1,i-1)).eq.0.0d0) then
        	uap=0.0d0
     	else
        	uap=(dsqrt(pv(1,i))*pv(2,i) + dsqrt(pv(1,i-1))*pv(2,i-1))/(dsqrt(pv(1,i))   + dsqrt(pv(1,i-1)))
     	endif
     	cap = dsqrt(g*0.5d0*(pv(1,i-1)+pv(1,i)))
     	helpmax=dmax1(dabs(uap-cap),dabs((uap+cap)))
     	if(helpmax.gt.maxl) then
       		 maxl=helpmax
     	endif
	mass=mass+pv(1,i)
ENDDO
dt = (cfl*dx)/maxl

!write(36,*) t,(engn_tot0-engn_tot)/engn_tot0,engn_tot,engn_tot0

me = (massin-mass)/massin
write(35,*) t,me


if(t.lt.Tt) then
	if(itert2.ge.0.02d0 ) then
       		call wsg(q,gp,t,be)
       		itert2=0.0d0
    	endif
    	itert2=itert2+dt

   	if(itert.ge.plt) then
   	!if(t.ge.pltt(pt) .and.pt.le.10) then
   	!if(t.ge.pltt(pt) .and.pt.le.13) then
   		it=it+1
   		call popr(q,be,dxi,g,t,it) !post processing
      		pt=pt+1
       		itert=0.0d0
   	endif
   	if(t+dt.gt.Tt) then
      		dt = Tt-t
   	endif
   	iter=iter+1
   	itert=itert+dt
   	goto 100
else
   	print*,'TELOS XRONOU '
endif
110  return
      end
