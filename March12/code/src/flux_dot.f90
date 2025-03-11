subroutine flx_dot(qp,pv,g,fl,entropy,order,limiter,sh,bed,bb,fh,nutot)
implicit none 
integer n,i,entropy,order,CL,CR,limiter,j,ii,k,l
character side,s
integer iwd(0:n+1),tmporcl,tmporcr,bb(n)
double precision pv(2,-3:n+4),bed(-3:n+4),swd(-3:n+4),eta(0:n+1)
!double precision qp(2,-3:n+4),flp(2,0:n+1),flm(2,0:n+1),f(2,0:n+1)
double precision qp(2,-3:n+4),fl(2,0:n+1),f(2,0:n+1)
double precision sh(2,n),duml(2,0:n),dumr(2,0:n)
double precision corr(0:n+1),corl(0:n+1),cor(0:n+1)
double precision uap,cap,lap1,lap2,a1,a2,rap1a
double precision rap1b,rap2a,rap2b,b1,b2
double precision w1a,w1b,w2a,w2b,sumq(2),g
double precision zm,up,um,u,tmp,c,cm,cp,hap,l1,l2
double precision f1abs,f2abs,tmpl,tmpr,dh,dhu
double precision uas11,uas12,uas13,uas21,uas22,uas23
double precision uas21p,uas22p,uas23p
double precision CCL,CCR,DB,threshold,SL,SR,lap11,lap22
double precision dt,dx,fh(0:n),nutot(n),energyf(n)
double precision R(2,2),R_inv(2,2),Cs(2,2), Ctotal(2,2)
double precision sgl(3),wgl(3),Integral,DQ(2),Jacob(2,2)
double precision ph1,ph2
common /tr/ threshold
common /dn/ dt,dx,n

    INTERFACE
        FUNCTION matmul_2x2(A, B) RESULT(C)
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: A(2,2), B(2,2)
            REAL(8) :: C(2,2)
        END FUNCTION matmul_2x2
        FUNCTION matvecmul_2x2(A, B) RESULT(C)
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: A(2,2), B(2)
            REAL(8) :: C(2)
        END FUNCTION matvecmul_2x2
    END INTERFACE

sgl(1) = (5.0d0 +dsqrt(15.0d0))/10.0d0
sgl(2) = 0.5d0
sgl(3) = (5.0d0 -dsqrt(15.0d0))/10.0d0
wgl(1) = 5.0/18.0d0
wgl(2) = 8.0d0/18.0d0
wgl(3) = wgl(1)

DO I=0,N  !faces 

        CL=I ! left cell
        CR=I+1! right cell
        uas11=pv(1,CL)   
        uas12=pv(2,CL) ! 1* L --- 2* R
        uas21=pv(1,CR)
        uas22=pv(2,CR)

        do k=1,2
             do l=1,2
                Ctotal(k,l) = 0.0d0
             enddo
        enddo
!--------------------------------------------------------------------------
        !Compute the Jacobial evaluated in the path for each sgl
        do ii=1,3
                ph1 = qp(1,cl) +sgl(ii)*(qp(1,cr)-qp(1,cl))
                ph2 = qp(2,cl) +sgl(ii)*(qp(2,cr)-qp(2,cl))

                uap = ph2/ph1
                hap = ph1

                cap = dsqrt(g*hap)
                lap1 = uap - cap
                lap2 = uap + cap

                Jacob(1,1) = 0.0d0
                Jacob(1,2) = 1.0d0
                Jacob(2,1) = cap**2-uap**2
                Jacob(2,2) = 2.0d0*uap

                R(1,1)= 1.0d0
                R(2,1)= lap1
                R(1,2)= 1.0d0
                R(2,2)= lap2

                !inverse matrix already multiplied with |\Lambda|
                R_inv(1,1)= dabs(lap1)*(1.0d0 +lap1/(2.0*cap))
                R_inv(2,1)= dabs(lap2)*(-lap1/(2.0d0*cap))
                R_inv(1,2)= dabs(lap1)*(-1.0d0/(2.0d0*cap))
                R_inv(2,2)= dabs(lap2)*(1.0d0/(2.0d0*cap))

                Cs(:,:)=0.0d0
                !Cs = matmul_2x2(R,R_inv)
                Cs = MATMUL(R,R_inv)

                do k=1,2
                   do l=1,2
                      Ctotal(k,l) = Ctotal(k,l) + wgl(ii)*Cs(k,l)!this is the integral
                   enddo
                enddo
        enddo        


        f(1,CL)= uas11*uas12
        f(2,CL)= (uas11*uas12**2) + 0.5d0*g*(uas11**2)
        f(1,CR)= uas21*uas22
        f(2,CR)= (uas21*uas22**2) + 0.5d0*g*(uas21**2)

        DQ(1) = qp(1,cr)-qp(1,cl)
        DQ(2) = qp(2,cr)-qp(2,cl)

        sumq = matvecmul_2x2(Ctotal,DQ)


        fl(1,i) =  0.5d0*(f(1,cl)+f(1,cr))-(0.50d0*sumq(1))
        fl(2,i) =  0.5d0*(f(2,cl)+f(2,cr))-(0.50d0*sumq(2))

ENDDO

DO I=1,N
        SH(1,I) = 0.0d0 
        SH(2,I) = 0.0d0 
ENDDO

return
end

