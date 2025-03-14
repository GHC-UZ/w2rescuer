subroutine flx_dot(qp,pv,g,flp,flm,entropy,order,limiter,sh,bed,bb,fh,nutot,por)
implicit none
integer n,i,entropy,order,CL,CR,limiter,j,ii,k,l
character side,s
integer iwd(0:n+1),tmporcl,tmporcr,bb(n)
double precision pv(4,-3:n+4),bed(-3:n+4),swd(-3:n+4),eta(-3:n+4)
double precision qp(4,-3:n+4),flp(4,0:n+1),flm(4,0:n+1),f(4,0:n+1)
double precision sh(4,n),por(-3:n+4)
double precision corr(0:n+1),corl(0:n+1),cor(0:n+1)
double precision uap,cap,lap1,lap2,lap3,lap4
double precision sumqp(4),sumqm(4),g
double precision zm,up,um,u,tmp,c,cm,cp,hap,l1,l2
double precision f1abs,f2abs,tmpl,tmpr,dh,dhu
double precision uas11,uas12,uas13,uas21,uas22,uas23
double precision etaap
double precision threshold
double precision dt,dx,fh(0:n),nutot(n),energyf(n)
double precision R(4,4),R_inv(4,4),Cs(4,4), Ctotalp(4,4),Ctotalm(4,4)
double precision sgl(4),wgl(4),Integral,DQ(4),Jacob(4,4)
double precision ph1,ph2,ph3,ph4,UL,UR
common /tr/ threshold
common /dn/ dt,dx,n

sgl(1) = (5.0d0 +dsqrt(15.0d0))/10.0d0
sgl(2) = 0.5d0
sgl(3) = (5.0d0 -dsqrt(15.0d0))/10.0d0
wgl(1) = 5.0/18.0d0
wgl(2) = 8.0d0/18.0d0
wgl(3) = wgl(1)

!sgl(1) = 0.5d0
!wgl(1) = 1.0d0

do i=-3,n+4
   eta(i) = qp(1,i)/qp(3,i)
enddo



DO I=0,N  !faces

        CL=I ! left cell
        CR=I+1! right cell

        UL = qp(2,cl)/(qp(1,cl)-qp(3,cl)*qp(4,cl)) 
        UR = qp(2,cr)/(qp(1,cr)-qp(3,cr)*qp(4,cr)) 

        do k=1,4
             do l=1,4
                Ctotalp(k,l) = 0.0d0
                Ctotalm(k,l) = 0.0d0
             enddo
        enddo
!--------------------------------------------------------------------------
        !Compute the Jacobian evaluated in the path for each sgl
        do ii=1,3
                ph1 = qp(1,cl) +sgl(ii)*(qp(1,cr)-qp(1,cl))
                ph2 = qp(2,cl) +sgl(ii)*(qp(2,cr)-qp(2,cl))
                ph3 = qp(3,cl) +sgl(ii)*(qp(3,cr)-qp(3,cl))
                ph4 = qp(4,cl) +sgl(ii)*(qp(4,cr)-qp(4,cl))

                uap = ph2/(ph1-ph3*ph4)
                etaap = ph1/ph3
                hap = etaap - ph4


                cap = dsqrt(g*hap)
                lap1 = uap - cap
                lap2 = 0.0d0
                lap3 = 0.0d0
                lap4 = uap + cap


                Jacob(1,1) = 0.0d0
                Jacob(1,2) = 1.0d0
                Jacob(1,3) = 0.0d0
                Jacob(1,4) = 0.0d0
                Jacob(2,1) = cap**2-uap**2
                Jacob(2,2) = 2.0d0*uap
                Jacob(2,3) = uap**2*ph4-g*etaap**2+g*etaap*ph4
                Jacob(2,4) = uap**2*ph3
                Jacob(3,1) = 0.0d0
                Jacob(3,2) = 0.0d0
                Jacob(3,3) = 0.0d0
                Jacob(3,4) = 0.0d0
                Jacob(4,1) = 0.0d0
                Jacob(4,2) = 0.0d0
                Jacob(4,3) = 0.0d0
                Jacob(4,4) = 0.0d0


                R(1,1)= 1.0d0
                R(2,1)= lap1
                R(3,1)= 0.0d0
                R(4,1)= 0.0d0
                R(1,2)= -Jacob(2,3)/Jacob(2,1)
                R(2,2)= 0.0d0
                R(3,2)= 1.0d0
                R(4,2)= 0.0d0 
                R(1,3)= -Jacob(2,4)/Jacob(2,1)
                R(2,3)= 0.0d0
                R(3,3)= 0.0d0
                R(4,3)= 1.0d0
                R(1,4)= 1.0d0
                R(2,4)= lap4
                R(3,4)= 0.0d0
                R(4,4)= 0.0d0

                !inverse matrix already multiplied with |\Lambda|
                R_inv(1,1)= dabs(lap1)*R(2,4)/(-R(2,1)+R(2,4))
                R_inv(2,1)= 0.0d0 
                R_inv(3,1)= 0.0d0
                R_inv(4,1)= dabs(lap4)*R(2,1)/(R(2,1)-R(2,4))
                R_inv(1,2)= dabs(lap1)*(-1.0d0)/(-R(2,1)+R(2,4))
                R_inv(2,2)= 0.0d0
                R_inv(3,2)= 0.0d0
                R_inv(4,2)= dabs(lap4)*1.0/(-R(2,1)+R(2,4))
                R_inv(1,3)= dabs(lap1)*R(1,2)*R(2,4)/(R(2,1)-R(2,4))
                R_inv(2,3)= 0.0d0*1.0d0
                R_inv(3,3)= 0.0d0
                R_inv(4,3)= dabs(lap4)*(-R(1,2)*R(2,1))/(R(2,1)-R(2,4))
                R_inv(1,4)= dabs(lap1)*(R(1,3)*R(2,4))/(R(2,1)-R(2,4))
                R_inv(2,4)= 0.0d0
                R_inv(3,4)= 0.0d0*1.0d0
                R_inv(4,4)= dabs(lap4)*(-R(1,3)*R(2,1))/(R(2,1)-R(2,4))

                Cs(:,:)=0.0d0
                Cs = MATMUL(R,R_inv)

                do k=1,4
                   do l=1,4
                      Ctotalp(k,l) = Ctotalp(k,l) + wgl(ii)*Cs(k,l) +wgl(ii)*Jacob(k,l)!this is the integral
                      Ctotalm(k,l) = Ctotalm(k,l) - wgl(ii)*Cs(k,l) +wgl(ii)*Jacob(k,l)!this is the integral
                   enddo
                enddo
        enddo !loop over ii



        DQ(1) = qp(1,cr)-qp(1,cl)
        DQ(2) = qp(2,cr)-qp(2,cl)
        DQ(3) = qp(3,cr)-qp(3,cl)
        DQ(4) = qp(4,cr)-qp(4,cl)


        sumqp = MATMUL(Ctotalp,DQ)
        sumqm = MATMUL(Ctotalm,DQ)


        flp(1,i) =  (0.50d0*sumqp(1))
        flp(2,i) =  (0.50d0*sumqp(2))
        flp(3,i) =  (0.50d0*sumqp(3))
        flp(4,i) =  (0.50d0*sumqp(4))

        flm(1,i) =  (0.50d0*sumqm(1))
        flm(2,i) =  (0.50d0*sumqm(2))
        flm(3,i) =  (0.50d0*sumqm(3))
        flm(4,i) =  (0.50d0*sumqm(4))
ENDDO !loop over the faces

return
end

