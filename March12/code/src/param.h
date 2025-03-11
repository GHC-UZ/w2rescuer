      implicit none
      integer n,dime,nmax,itime
      parameter(nmax=80000,dime=2)
      double precision g,pow4,sx,tmpx,theta,Nm
      double precision x,dx,dt,a,b,cfl,t,Tt,dt1,dt2
      double precision q(2,-3:nmax+4),qp(2,-3:nmax+4),q1(2,-3:nmax+4)
      double precision pv(2,-3:nmax+4),pr(2,-3:nmax+4),SWD(-3:NMAX+4)
      double precision qs(2,-3:nmax+4),qphi(2,-3:nmax+4)
      double precision PHI1(NMAX),PHI2(NMAX),PHIP(NMAX)
      double precision F(2,NMAX),F1(2,NMAX),F2(2,NMAX),FP(2,NMAX)
      double precision QM(2,-3:NMAX+4),FT(2,NMAX)
      double precision K1,K2
      double precision DL(NMAX),D(NMAX),DU(NMAX)
      double precision fs(dime,nmax),maxl,fhm(2,nmax),fhs(2,nmax)
      double precision fm(dime,nmax),fl(dime,0:nmax+1)
      double precision flp(dime,nmax),flm(dime,0:nmax+1)
      double precision helpmax,cap,tmp,maxv
      double precision sh(2,nmax),vv(8,nmax)
      double precision s(2,nmax),be(-3:nmax+4),itert2
      double precision u,up,uap,itert
      integer i,entropy,method,limiter,hr,problem,iter,order,gp(600)
      integer flag,flagtmp,case,tm,it,ITE,INFO,IRK,pt,flagm,iwagen
      double precision mass,me,massin,etot,threshold,hd,plt,rb(nmax)
      double precision dxi(0:nmax+1),xc(0:nmax+1),df,error,fr(nmax)
      double precision qstar(2,-3:nmax+4),lsl,lsr,hh,error2
      double precision t1,t2,t3,wagen,swl
      double precision tmpt,fh(0:nmax),Nutot(nmax),qsponge(2,nmax)
      double precision PHI(nmax),PHIB(nmax),W(nmax),R(nmax),abar
      double precision engn_tot0,engn_tot,energyf(nmax),amplitude
      integer BB(nmax),choisebr,fdt,iflux
      common /sc/ cfl,Tt,plt,lsl,lsr,order,entropy,limiter,fdt
      common /da/hd,case
      common /dn/ dt,dx,n
      common /bd/ a,b
      common /tr/ threshold
      common /fr/ Nm,theta
      common /gen/ tmpt,iwagen
      common /gnpar/ abar,choisebr
      common /energy/engn_tot0,engn_tot,energyf 
      common /flux/ iflux

