subroutine friction(Q,PV,g,fr)
implicit none
integer n,i
double precision dt,dx,Nm,theta,threshold,g
double precision Q(2,-3:N+4),PV(2,-3:N+4),fr(n)
double precision RFP,RF,vel,dum1,dum2,up,tmp
common /dn/ dt,dx,n
common /tr/ threshold
common /fr/ Nm,theta

!---------------EXPLICIT TREATMENT---------------
!do i=1,n
! 	if(pv(1,i).le.1.d-4 .or. Nm.eq.0.0d0) then
!  	!if(pv(1,i).le.threshold .or. Nm.eq.0.0d0) then
!   		DUM1=0.0D0
!  	else
!   		RFP=Nm*Nm*dabs(pv(2,i))/(q(1,i)*dexp(dlog(q(1,i))/3.0d0))
!   		DUM1 = (-g*pv(2,i)*q(1,i)*RFP)
!  	endif
!	fr(i)=DUM1
!enddo
!
!-------------------semi implicit-----------------
do i=1,n
  	if(pv(1,i).le.1.d-4) then
   		DUM1=0.0D0
 	else
   		RFP=Nm*Nm*dabs(pv(2,i))/(pv(1,i)*dexp(dlog(pv(1,i))/3.0d0))
   		DUM1 = (-g*theta*pv(2,i)*pv(1,i)*RFP*dt)
  	endif

  	if(q(1,i).le.1.d-4) then
   		DUM2=1.0d0
  	else
   		vel=q(2,i)/q(1,i)
   		RF =Nm*Nm*dabs(vel)/(dexp(dlog(q(1,i))/3.0d0)*q(1,i))
   		DUM2 = (1.0d0+(1.0d0-theta)*g*RF*dt)
  	endif
   	Q(2,I) = (Q(2,I)+DUM1)/DUM2
enddo
return
end



