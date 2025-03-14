subroutine initial(qp,n,dx,a,g,maxl,sum,etot,b,dxi,amplitude,swl,por)
implicit none
integer n,problem,i1,i2,i3,i,nh,k,is,ie,case,hr,entropy,limiter
integer ic,order
double precision qp(4,-3:n+4),dx,a,maxl,helpmax,ls,por(-3:n+4)
double precision g,B(-3:n+4),ho,A1,swd(-3:n+4),pv(4,-3:n+4)
double precision bin(-3:n+4),bfn(-3:n+4),t
double precision dt,sum,u,up,uap,cap,swl
double precision x,xd(1001),hta(1001),hta1(n),tmp,tmp1,etot
double precision dxi(0:n+1),ircount
double precision gamma,x1,pi,hd,cfl,Tt,plt,amplitude
common /sc/ cfl,Tt,plt,ls,order,entropy,limiter
common /da/hd,case

t=0.0d0
x=a
maxl=0.0d0
tmp=dx
tmp1=dx!/2.0d0
sum=0.0d0
!-------------------------------------------------------------------------
!--------------- CONSTRACTION OF THE INITIAL BED--------------------------
!-------------------------------------------------------------------------
call pr1_bed(b,dxi)
!-------------------------------------------------------------------------
!----------------- INITIAL CONDITIONS FOR H AND U -------------------------
!-------------------------------------------------------------------------
call pr1_initial(qp,b,0.0d0,amplitude,swl,por)
!************************************************************************* 
call primitives(qp,pv,por,b)
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
ENDDO

do i=1,n-1
         sum=sum+qp(1,i)*dx*0.5d0+qp(1,i+1)*dx*0.5d0
enddo

return
end

