subroutine pr1_bed(be,dxi)
implicit none
integer n,i, dias, regular
double precision be(-3:n+4), c1, c2, y0, y1, x0, x1
double precision x,dt,dx,dxi(0:n+1),h_an,a,b
common /dn/ dt,dx,n
common /bd/ a,b

!--------------------------CONSTANT------------------------------
do i=1,n
       be(i)=0.0d0
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
      
