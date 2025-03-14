subroutine pr1_bed(be,dxi)
implicit none
integer n,i, dias, regular
double precision be(-3:n+4), c1, c2, y0, y1, x0, x1
double precision x,dt,dx,dxi(0:n+1),h_an,a,b,pi
common /dn/ dt,dx,n
common /bd/ a,b

pi=4.0d0*datan(1.0d0)
!--------------------------CONSTANT------------------------------
x=a
do i=1,n
       be(i)=0.0d0
       !be(i) =  0.05d0*dsin(x-12.5)*dexp(1.0d0-(x-25.d0)*(x-25.d0))
       !if(1.4d0 <= x .and. x<= 1.6d0) then
       !        be(i) = -0.25d0*(dcos(10.0d0*pi*(x-15.0d0)) -1.0d0)
       !else
       !        be(i)=0.0d0
       !endif
       x=x+dx
enddo

call bdc_bed(be)

return
end
!------------------------------------------------------------------------------
subroutine bdc_bed(b) 
implicit none
integer n,problem
double precision b(-3:n+4),dt,dx
common /dn/dt,dx,n

  B(0) = B(2)
  B(-1) = B(3)
  B(-2) = B(4)
  B(-3) = B(5)
  B(n+1) = B(n-1)
  B(n+2) = B(n-2)
  B(n+3) = B(n-3)
  B(n+4) = B(n-4)

return
end
      
