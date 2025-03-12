subroutine euler(qp,pv,be,gp,ite,massin,por)
include 'param.h'
double precision tmp_phi(0:n+1),PLTT(10)


if(Tt.eq.0.0d0) goto 110
g=9.81d0
iter=1
itert=0.0d0
ITE=0
it=0
100  t = t + dt     

flagtmp=1
mass=0.0d0
etot=0.0d0

      print*,'TIME: ' ,t , ' dt: ',dt    
      ITE=ITE+1
!*************************************************************************

x=a
flag=0
maxl=0.0d0
helpmax=0.0d0
maxv=0.0d0

!-----everything is computed using values from the previous time step
call primitives(qp,pv,por,be)
call flx_dot(qp,pv,g,flp,flm,entropy,order,limiter,sh,be,bb,fh,nutot,por) !calculate the fluxes ADVECTIVE PART

!############################## Implement here the fluxes ############################
!---------------------------------------FIND H----------------------------------
DO I=1,N
        Q(1,I) = QP(1,I) + DT*F(1,I) 
        Q(2,I) = QP(2,I) + DT*F(2,I) 
        Q(3,I) = QP(3,I) + DT*F(3,I) 
        Q(4,I) = QP(4,I) + DT*F(4,I) 
ENDDO

call friction(q,pv,g,fr)

call boundarya(q,n,t)
call boundaryb(q,n)

DO I=-3,N+4
        QP(1,I)=Q(1,I)
        QP(2,I)=Q(2,I)
        QP(3,I)=Q(3,I)
        QP(4,I)=Q(4,I)
ENDDO 

!      call runup(q,be,n,t)
!************************************************************************
call primitives(q,pv,por,be)
mass=0.0d0
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
      	if(maxv.lt.dabs(vv(1,i)).and. rb(i).ne.0.0d0) maxv=dabs(vv(1,i))
	mass=mass+pv(1,i)
ENDDO
dt = (cfl*dx)/maxl


me = (massin-mass)/massin
write(35,*) t,me
if(t.lt.Tt) then
	 if(itert2.ge.0.02d0 ) then
    		call wsg(q,gp,t,be)
    		itert2=0.0d0
 	endif
 	itert2=itert2+dt

	if(itert.ge.plt) then
		it=it+1
		call popr(q,be,dxi,g,t,it,por) !post processing
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
