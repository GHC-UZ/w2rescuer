subroutine gauge(gp)
implicit none
integer i,j,n,gp(600),nm
double precision x
double precision dt,dx,a,b,g(600),dis(600)
common/dn/ dt,dx,n
common/bd/ a,b

nm=0
open(31,file='../input/gauges.dat')
do i=1,599
	nm=nm+1
   	read(31,*,end=100) g(i)
   	gp(i)=0
	dis(i)=1000.0d0
enddo
100  close(31)
gp(600)=nm-1

x=a
do j=1,n
	do i=1,nm-1
      		if(dabs(x-g(i)).lt.dis(i)) then  
         		gp(i)= j
         		dis(i)=dabs(x-g(i))
      		endif
   	enddo
   	x=x+dx
enddo
return
end
!***********************************************************************
subroutine wsg(q,gp,t,be)
implicit none
integer i,gp(600),n,nm
double precision be(-3:n+4)
double precision dt,dx,q(2,-3:n+4),t,swd(-3:n+4)
common/dn/ dt,dx,n

nm=gp(600)
!write(33,*) t,(q(1,gp(i))-swd(gp(i)), i=1,nm)
write(32,*) t,(q(1,gp(i))+be(gp(i)), i=1,nm)

return
end
