      subroutine muscl(uas11,uas12,uas13,uas21,uas22,uas23,pv,bed,cl,cr,order,lim,iwd,bb,fh,nutot)
      implicit none
      integer n,i,order,CL,CR,lim,l
      integer iwd(0:n+1),bb(n)
      double precision pv(2,-3:n+4),bed(-3:n+4),threshold
      double precision uas11,uas12,uas13,uas21,uas22,uas23
      double precision fhi,dh,du,dhm,dhp,dum,dup,r1,r2,r11,r22
      double precision k,dsh,dsu,dshm,dshp,dsum,dsup,eps,dx,dt
      double precision db,dsb,dbm,dbp,dsbm,dsbp,dhav,duav,dbav
      double precision fh(0:n),nutot(n),energyf
      common /dn/ dt,dx,n
      common /tr/threshold

eps=1.d-15
if(order.eq.3) then 
   	k=1.0d0/3.0d0
elseif(order.eq.2) then 
   	k=-1.0d0 !fully upwind scheme
  	!k=0.0d0 !fromm scheme
  	!k=1.0d0 !three-point central-difference scheme 
endif

l=lim
!-------LEFT---------------------
dhm=dh(pv,CL-1,CR-1,n) !dh_{i-1/2}
dhp=dh(pv,CL,CR,n) !dh_{i+1/2}
dum=du(pv,CL-1,CR-1,n)
dup=du(pv,CL,CR,n)
dbm=db(bed,CL-1,CR-1,n)
dbp=db(bed,CL,CR,n)

!---------------------2order------------------------------
if(order.eq.2) then 
	k=1.0d0
   	r1 = dhm/(dhp+eps)
   	uas11 =pv(1,CL)+fhi(fh,cl,cr,r1,l,bb)*((1.0d0-k)*dhm+(1.0d0+k)*dhp)/4.0d0
  	r1 = dum/(dup+eps)
   	uas12 =pv(2,CL)+fhi(fh,cl,cr,r1,l,bb)*((1.0d0-k)*dum+(1.0d0+k)*dup)/4.0d0
   	r1 = dbm/(dbp+eps)
   	uas13 =bed(CL)+fhi(fh,cl,cr,r1,l,bb)*((1.0d0-k)*dbm+(1.0d0+k)*dbp)/4.0d0

  	if(pv(1,CL).GT.0.0D0 .and. pv(1,CL-1).EQ.0.0D0 ) IWD(CL)=1
  	if(pv(1,CL).GT.0.0D0 .and. pv(1,CL-2).EQ.0.0D0 ) IWD(CL)=1
  	if(pv(1,CL).GT.0.0D0 .and. pv(1,CR).EQ.0.0D0 ) IWD(CL)=1
endif
!---------------------3 order------------------------------
if(order.eq.3) then 
	r1 = dhm/(dhp+eps)
    	r11=1.0d0/(r1+eps)
    	if(r1.eq.0.0d0) r11=0.0d0
    	uas11 = pv(1,CL) +( (1.0d0-k)*dhm*fhi(fh,cl,cr,r1,l,bb) +(1.0d0+k)*dhp*fhi(fh,cl,cr,r11,l,bb) )/4.0d0
    	r1 = dum/(dup+eps)
    	uas12 = pv(2,CL) +( (1.0d0-k)*dum*fhi(fh,cl,cr,r1,l,bb) +(1.0d0+k)*dup*fhi(fh,cl,cr,r11,l,bb) )/4.0d0
    	r1 = dbm/(dbp+eps)
    	uas13 = bed(CL) +( (1.0d0-k)*dbm*fhi(fh,cl,cr,r1,l,bb) +(1.0d0+k)*dbp*fhi(fh,cl,cr,r11,l,bb) )/4.0d0

!---check if a dry node is involved
        if(pv(1,CL).GT.0.0D0 .and. pv(1,CL-1).EQ.0.0D0 ) IWD(CL)=1
        if(pv(1,CL).GT.0.0D0 .and. pv(1,CL-2).EQ.0.0D0 ) IWD(CL)=1
        if(pv(1,CL).GT.0.0D0 .and. pv(1,CR).EQ.0.0D0 ) IWD(CL)=1
endif

!-------RIGHT------------------
      dhm=dh(pv,CL,CR,n)  !=previous dhp d_{i+1/2}
      dhp=dh(pv,CL+1,CR+1,n) ! dh_{i+3/2}
      dum=du(pv,CL,CR,n)
      dup=du(pv,CL+1,CR+1,n)
      dbm=db(bed,CL,CR,n)
      dbp=db(bed,CL+1,CR+1,n)
!
!---------------------2 order------------------------------
if(order.eq.2) then 
     k=-1.0d0
     r2 = dhm/(dhp+eps)
     uas21 = pv(1,CR)-fhi(fh,cl,cr,r2,l,bb)*((1.0d0-k)*dhp+(1.0d0+k)*dhm)/4.0d0
     r2 = dum/(dup+eps)
     uas22 = pv(2,CR)-fhi(fh,cl,cr,r2,l,bb)*((1.0d0-k)*dup+(1.0d0+k)*dum)/4.0d0
     r2 = dbm/(dbp+eps)
     uas23 = bed(CR)-fhi(fh,cl,cr,r2,l,bb)*((1.0d0-k)*dbp+(1.0d0+k)*dbm)/4.0d0

     if(pv(1,CR).GT.0.0D0 .and. pv(1,CR-1).EQ.0.0D0 ) IWD(CR)=1
     if(pv(1,CR).GT.0.0D0 .and. pv(1,CR+1).EQ.0.0D0 ) IWD(CR)=1
     if(pv(1,CR).GT.0.0D0 .and. pv(1,CR+2).EQ.0.0D0 ) IWD(CR)=1
  endif

!---------------------3 order------------------------------
if(order.eq.3) then 
	r2 = dhm/(dhp+eps)
 	r22=1.0d0/(r2+eps)
 	if(r2.eq.0.0d0) r22=0.0d0
	uas21 = pv(1,CR) -((1.0d0-k)*dhp*fhi(fh,cl,cr,r2,l,bb)+(1.0d0+k)*dhm*fhi(fh,cl,cr,r22,l,bb))/4.0d0

	r2 = dum/(dup+eps)
	r22=1.0d0/r2
	if(r2.eq.0.0d0) r22=0.0d0
 	uas22 = pv(2,CR) -((1.0d0-k)*dup*fhi(fh,cl,cr,r2,l,bb)+(1.0d0+k)*dum*fhi(fh,cl,cr,r22,l,bb))/4.0d0

	r2 = dbm/(dbp+eps)
	r22=1.0d0/r2
 	if(r2.eq.0.0d0) r22=0.0d0
 	uas23 = bed(CR) -((1.0d0-k)*dbp*fhi(fh,cl,cr,r2,l,bb) +(1.0d0+k)*dbm*fhi(fh,cl,cr,r22,l,bb))/4.0d0

 	if(pv(1,CR).GT.0.0D0 .and. pv(1,CR-1).EQ.0.0D0 ) IWD(CR)=1
 	if(pv(1,CR).GT.0.0D0 .and. pv(1,CR+1).EQ.0.0D0 ) IWD(CR)=1
 	if(pv(1,CR).GT.0.0D0 .and. pv(1,CR+2).EQ.0.0D0 ) IWD(CR)=1
endif


if(uas11.LE.threshold) then 
   	uas11=pv(1,CL)
   	uas12=pv(2,CL)
   	uas13=bed(CL)
endif
if(uas21.LE.threshold) then 
   	uas21=pv(1,CR)
   	uas22=pv(2,CR)
   	uas23=bed(CR)
endif

return 
end
!------------------------------------------------------------------------------
double precision function dh(pv,CL,CR,N)
implicit none
integer N,CL,CR
double precision pv(2,-3:n+4)
 dh = pv(1,CR)-pv(1,CL)
return 
end
!------------------------------------------------------------------------------
double precision function du(pv,CL,CR,N)
implicit none
integer CL,CR,N
double precision pv(2,-3:n+4)
 du = pv(2,CR) -pv(2,CL)
return 
end
!------------------------------------------------------------------------------
double precision function db(bed,CL,CR,N)
implicit none
integer CL,CR,N
double precision bed(-3:n+4)
 db = bed(CR) -bed(CL)
return 
end
!------------------------------------------------------------------------------
double precision function dsh(pv,CL,CR,N)
implicit none
integer N,CL,CR
double precision pv(2,-3:n+4),minmod3
double precision dum1,dum2,dum3,d3,dh
double precision dum11,dum22,dum33

dum1=dh(pv,CL,CR,N) !dh_{i+1/2} cl=i
dum2=dh(pv,CL+1,CR+1,N) !dh_{i+3/2}
dum3=dh(pv,CL-1,CR-1,N) !dh_{i-1/2}

dum11=minmod3(dum1,dum2,dum3)
dum22=minmod3(dum2,dum3,dum1)
dum33=minmod3(dum3,dum1,dum2)

d3 = dum22-2.0d0*dum11+dum33
!d3 = dum2-2.0d0*dum1+dum3
dsh = dum1 -d3/6.0d0
return 
end
!-------------------------------------------------------------------------
double precision function dsu(pv,CL,CR,N)
implicit none
integer N,CL,CR
double precision pv(2,-3:n+4),minmod3
double precision dum1,dum2,dum3,d3,du
double precision dum11,dum22,dum33

dum1=du(pv,CL,CR,N) !dh_{i+1/2} cl=i
dum2=du(pv,CL+1,CR+1,N) !dh_{i+3/2}
dum3=du(pv,CL-1,CR-1,N) !dh_{i-1/2}

dum11=minmod3(dum1,dum2,dum3)
dum22=minmod3(dum2,dum3,dum1)
dum33=minmod3(dum3,dum1,dum2)

d3 = dum22-2.0d0*dum11+dum33
!d3 = dum2-2.0d0*dum1+dum3
dsu = dum1 -d3/6.0d0
return 
end
!--------------------------------------------------------------------------
double precision function dsb(bed,CL,CR,N)
implicit none
integer N,CL,CR
double precision bed(-3:n+4),minmod3
double precision dum1,dum2,dum3,d3,db
double precision dum11,dum22,dum33

dum1=db(bed,CL,CR,N) !dh_{i+1/2} cl=i
dum2=db(bed,CL+1,CR+1,N) !dh_{i+3/2}
dum3=db(bed,CL-1,CR-1,N) !dh_{i-1/2}

dum11=minmod3(dum1,dum2,dum3)
dum22=minmod3(dum2,dum3,dum1)
dum33=minmod3(dum3,dum1,dum2)

d3 = dum22-2.0d0*dum11+dum33
 !d3 = dum2-2.0d0*dum1+dum3
dsb = dum1 -d3/6.0d0
return 
end
!--------------------------------------------------------------------------
double precision function dhav(pv,CL,CR,N,l)
implicit none
integer N,CL,CR,l,bb(n)
double  precision pv(2,-3:n+4),eps
double precision dh,dum1,dum2,fhi,r,fh(0:n)

fh(:)=1.0d0
eps=1.d-15
bb(:)=0
dum1=dh(pv,CL-1,CR-1,N)
dum2=dh(pv,CL,CR,N)

r=dum1/(dum2+eps)
dhav = 0.5d0*fhi(fh,cl,cr,r,l,bb)*(dum1+dum2)

return
end
!--------------------------------------------------------------------------
double precision function duav(pv,CL,CR,N,l)
implicit none
integer N,CL,CR,l,bb(n)
double  precision pv(2,-3:n+4),eps
double precision du,dum1,dum2,fhi,r,fh(0:n)

fh(:)=1.0d0
eps=1.d-15
bb(:)=0
dum1=du(pv,CL-1,CR-1,N)
dum2=du(pv,CL,CR,N)

r=dum1/(dum2+eps)
duav = 0.5d0*fhi(fh,cl,cr,r,l,bb)*(dum1+dum2)

return
end
!--------------------------------------------------------------------------
double precision function dbav(bed,CL,CR,N,l)
implicit none
integer N,CL,CR,l,bb(n)
double  precision bed(-3:n+4),eps
double precision db,dum1,dum2,fhi,r,fh(0:n)

fh(:)=1.0d0
eps=1.d-15
bb(:)=0
dum1=db(bed,CL-1,CR-1,N)
dum2=db(bed,CL,CR,N)

r=dum1/(dum2+eps)
dbav = 0.5d0*fhi(fh,cl,cr,r,l,bb)*(dum1+dum2)

return
end
