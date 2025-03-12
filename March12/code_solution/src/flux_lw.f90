      subroutine flx_lw(qp,pv,g,fl,entropy,order,limiter,sh,bed,bb,fh,nutot)
      implicit none 
      integer n,i,entropy,order,CL,CR,limiter,j
      character side,s
      integer iwd(0:n+1),tmporcl,tmporcr,bb(n)
      double precision pv(2,-3:n+4),bed(-3:n+4)
      double precision qp(2,-3:n+4),fl(2,0:n+1),f(2,0:n+1)
      double precision sh(2,n),duml(2,0:n),dumr(2,0:n)
      double precision corr(0:n+1),corl(0:n+1),cor(0:n+1)
      double precision uap,cap
      double precision sum1,sum2,g
      double precision zm,up,um,u,tmp,c,cm,cp,hap,l1,l2
      double precision f1abs,f2abs,tmpl,tmpr,dh,dhu
      double precision uas11,uas12,uas13,uas21,uas22,uas23
      double precision uas21p,uas22p,uas23p
      double precision DB,threshold,a1,a2,b1,b2
      double precision dt,dx,fh(0:n),nutot(n)
      common /tr/ threshold
      common /dn/ dt,dx,n

      DO J=1,N
         IWD(J)=0
      ENDDO
      DO I=0,N  !faces 

        CL=I ! left cell
        CR=I+1! right cell
        uas11=pv(1,CL)   
        uas12=pv(2,CL) ! 1* L --- 2* R
        uas21=pv(1,CR)
        uas22=pv(2,CR)

        if(dsqrt(uas11) + dsqrt(uas21) .eq. 0.0d0) then
            uap=0.0d0
        else
            uap = 0.5d0*(uas12+uas22) 
        endif
        ! COMPUTE A^2 =Jacobian^2
        ! We choose A_{1+1/2}=A((Q_L+Q_R)/2)

        hap = 0.5d0*(uas11 + uas21)
        cap = dsqrt(g*hap)
        a1=cap**2-uap**2
        a2=2.0d0*uap
        b1=2.0d0*uap*a1
        b2=a1+4.0d0*uap**2

        dh=uas21-uas11
        dhu=uas21*uas22-uas11*uas12
       
        f(1,CL)= uas11*uas12
        f(2,CL)= (uas11*uas12**2) + 0.5d0*g*(uas11**2)
        f(1,CR)= uas21*uas22
        f(2,CR)= (uas21*uas22**2) + 0.5d0*g*(uas21**2)

        sum1 = a1*dh + a2*dhu
        sum2 = b1*dh + b2*dhu

        fl(1,i) = (0.50d0*(f(1,CL)+f(1,CR))) -(dt/dx)*(0.50d0*sum1)
        fl(2,i) = (0.50d0*(f(2,CL)+f(2,CR))) -(dt/dx)*(0.50d0*sum2)

      ENDDO

      DO I=1,N
         SH(1,I) = 0.0d0 
         SH(2,I) = 0.0d0 
      ENDDO

      return
      end

