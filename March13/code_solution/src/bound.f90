subroutine boundarya(qp,n,t)
integer n
double precision qp(4,-3:n+4),t

qp(1,1)=qp(1,2)
qp(1,0) = qp(1,2)
qp(1,-1) = qp(1,3)
qp(1,-2) = qp(1,4)
qp(1,-3) = qp(1,5)

qp(2,0) = qp(2,2)
qp(2,-1) = qp(2,3)
qp(2,-2) = qp(2,4)
qp(2,-3) = qp(2,5)

qp(3,0) = qp(3,2)
qp(3,-1) = qp(3,3)
qp(3,-2) = qp(3,4)
qp(3,-3) = qp(3,5)

qp(4,0) = qp(4,2)
qp(4,-1) = qp(4,3)
qp(4,-2) = qp(4,4)
qp(4,-3) = qp(4,5)

return
end
!-----------------------------------------------------------
subroutine boundaryb(qp,n)
integer n
double precision qp(4,-3:n+4)

qp(1,n+1) = qp(1,n-1)
qp(1,n+2) = qp(1,n-2)
qp(1,n+3) = qp(1,n-3)
qp(1,n+4) = qp(1,n-4)

qp(2,n+1) = qp(2,n-1)
qp(2,n+2) = qp(2,n-2)
qp(2,n+3) = qp(2,n-3)
qp(2,n+4) = qp(2,n-4)

qp(3,n+1) = qp(3,n-1)
qp(3,n+2) = qp(3,n-2)
qp(3,n+3) = qp(3,n-3)
qp(3,n+4) = qp(3,n-4)

qp(3,n+1) = qp(3,n-1)
qp(3,n+2) = qp(3,n-2)
qp(3,n+3) = qp(3,n-3)
qp(3,n+4) = qp(3,n-4)

 return
 end
