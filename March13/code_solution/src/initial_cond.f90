subroutine pr1_initial(qp,be,ttmp,amplitude,swl,por)
implicit none
integer n,i
double precision qp(4,-3:n+4),be(-3:n+4),dt,dx,swd(-3:n+4)
double precision qp2(4,-3:n+4),swl,por(-3:n+4)
double precision g,x,d,t,a,b,pv(4,-3:n+4),ttmp
double precision solh,solu,amplitude,x0, pi
double precision gamma
common /dn/ dt,dx,n
common /bd/ a,b

pi=4.0d0*datan(1.0d0)
g=9.81d0
x=a 
t=0.0d0 !Initial time
d=swl! The still water level water depth
x0 = 10.0 !location of the dam
      
do i=1,n
!---------------DAM BREAK---------------------------------
        if(x<= 50) then   !INITIAL CONDITION AT THE LEFT
                por(i) = 0.9d0
                qp(1,i) = por(i)*8.0d0 !(porosity*eta)
                qp(2,i)= 0.0d0! (porosity* u)
        else
                por(i) = 0.7d0    !INITIAL CONDITION AT THE RIGHT
                qp(1,i)=por(i)*3.0d0
                qp(2,i)= 0.0d0
        endif
        qp(3,i) = por(i)
        qp(4,i) = be(i)

!---------------for dry areas-----------------------           
         if(qp(1,i).le.0.0d0) then 
                qp(1,i)=0.0d0
                qp(2,i)=0.0d0
         endif
         x=x+dx
enddo

call boundarya(qp,n,t)
call boundaryb(qp,n)

por(-3)=por(4)
por(-2)=por(3)
por(-1)=por(2)
por(0)=por(1)

por(n+1) = por(n)
por(n+2) = por(n-1)
por(n+3) = por(n-2)
por(n+4) = por(n-3)

call primitives(qp,pv,por,be)

return
end
!------------------------------------------------------------------------------
double precision function solh(x,x0,amplitude,d0,g)
implicit none
integer i,n
double precision dt,dx,d0,x,x0,g
double precision d1,d2,celerity,k
double precision amplitude,sech,xx
common /dn/ dt,dx,n

d1=d0
d2=amplitude+d1
celerity=sqrt(g*d2)

xx=0.5d0*(x-x0)*sqrt(3.0d0*amplitude/(d2*d1**2))
! k=dsqrt(3.0d0*amplitude)/(2.0d0*d1*sqrt(d2))
!xx=k*(x-x0)

solh= d0+amplitude*sech(xx)**2

!if(x.le.0.0d0 .or. x.gt.100.0d0) solh=d0
return
end
!------------------------------------------------------------------------------
double precision function solu(x,x0,amplitude,d0,g)
implicit none
integer i,n
double precision dt,dx,d0,x,x0,g
double precision d1,d2,celerity
double precision amplitude,solh
common /dn/ dt,dx,n

d1=d0
d2=amplitude+d1
celerity=sqrt(g*d2)

solu= celerity*(1.0d0-d1/solh(x,x0,amplitude,d0,g))

 !if(x.le.0.0d0 .or. x.gt.100.0d0) solu=0.0
return
end
!------------------------------------------------------------------------------
double precision function dsech(x)
implicit none
double precision sech,x,b

b=1.0d0/15.0d0
dsech = -b*sech(b*x)*dtanh(b*x)

return
end
!------------------------------------------------------------------------------
double precision function sech(x)
implicit none
double precision x

sech=1.0d0/dcosh(x)

return
end
!------------------------------------------------------------------------------
subroutine bc_ini(qp,b)
implicit none
integer n
double precision qp(2,-3:n+4),b(-3:n+4),dt,dx
common /dn/ dt,dx,n

qp(1,0)=qp(1,1)
qp(1,-1)=qp(1,2)
qp(1,-2)=qp(1,3)
qp(1,-3)=qp(1,4)
qp(2,0)=qp(2,1)
qp(2,-1)=qp(2,2)
qp(2,-2)=qp(2,3)
qp(2,-3)=qp(2,4)

qp(1,n+1) = qp(1,n)
qp(1,n+2) = qp(1,n-1)
qp(1,n+3) = qp(1,n-2)
qp(1,n+4) = qp(1,n-3)
qp(2,n+1) = qp(2,n)
qp(2,n+2) = qp(2,n-1)
qp(2,n+3) = qp(2,n-2)
qp(2,n+4) = qp(2,n-3)

return
end
