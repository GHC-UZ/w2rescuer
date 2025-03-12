subroutine grid(xi,problem,xc)
implicit none
integer n,problem,i,j,sn,en,midle,m
double precision a,b,dt,dx,xi(0:n+1),x,dx1
double precision sum,tmp(n/2),p
double precision xc(0:n+1),la,lb
common /dn/dt,dx,n
common /bd/a,b

x=a
do i=1,n
	xi(i) = x
   	x=x+dx
enddo         

xi(0) = xi(1)-dx
xi(n+1) = xi(n)+dx

do i=1,n
  	xc(i) = 0.5d0*(xi(i+1)+xi(i))
enddo
xc(0) = xc(1)-dx
xc(n+1) = xc(n)+dx


return
end
