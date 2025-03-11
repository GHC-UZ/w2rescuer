      subroutine flx_hll(qp,pv,g,fl,entropy,order,limiter,sh,bed,bb,fh,nutot)
      implicit none 
      integer n,i,entropy,order,CL,CR,limiter
      integer iwd(0:n+1),bb(n)
      character side,s
      double precision pv(2,-3:n+4),threshold
      double precision qp(2,-3:n+4),fl(2,0:n+1),f(2,0:n+1)
      double precision uas11,uas12,uas21,uas22,g,sh(2,n),bed(-3:n+4)
      double precision CCL,CCR,SL,SR,us,fs,HST,UST,QL,QR,SP
      double precision uas13,uas23,hta(2,-3:n+4)
      double precision uap,cap,hap,lap1,lap2,l1,l2,eps
      double precision dt,dx,fh(0:n),nutot(n),energyf(n)
      common /tr/ threshold
      common /dn/ dt,dx,n

      eps=1.d-15

      DO I=0,N  !faces 

        CL=I ! left cell
        CR=I+1! right cell
        uas11=pv(1,CL)   
        uas12=pv(2,CL) ! 1* L --- 2* R
        uas21=pv(1,CR)
        uas22=pv(2,CR)


        f(1,CL)= uas11*uas12
        f(2,CL)= (uas11*uas12**2) + 0.5d0*g*(uas11**2)
        f(1,CR)= uas21*uas22
        f(2,CR)= (uas21*uas22**2) + 0.5d0*g*(uas21**2)

        CCL=dsqrt(g*uas11)  
        CCR=dsqrt(g*uas21) 
        us=(uas12+uas22)/2.0d0 +CCL-CCR
        fs=(CCL+CCR)/2.0d0 +(uas12-uas22)/4.0d0
        SL=dmin1(uas12-CCL,us-fs)
        SR=dmax1(uas22+CCR,us+fs)


       IF(SL.GE.0.0D0) THEN
         fl(1,i)=f(1,CL)
         fl(2,i)=f(2,CL)
       ELSEIF(SL.LE.0.0D0 .AND. SR.GE.0.0D0)THEN
        fl(1,i) = (SR*f(1,CL)-SL*f(1,CR) +SL*SR*(uas21-uas11))/(SR-SL)
        fl(2,i) = (SR*f(2,CL)-SL*f(2,CR) +SL*SR*(uas22*uas21-uas12*uas11))/(SR-SL)
       ELSEIF(SR.LE.0.0D0) then
         fl(1,i)=f(1,CR)
         fl(2,i)=f(2,CR)
       ENDIF

      ENDDO

      DO I=1,N
         SH(1,I)=0.0D0
         SH(2,I)=0.0D0
      ENDDO

      return
      end
