subroutine flx_roe(qp,pv,g,fl,entropy,order,limiter,sh,bed,bb,fh,nutot)
implicit none 
integer n,i,entropy,order,CL,CR,limiter,j
integer slope
character side,s
integer iwd(0:n+1),tmporcl,tmporcr,bb(n)
double precision pv(2,-3:n+4),bed(-3:n+4),swd(-3:n+4),eta(0:n+1)
double precision qp(2,-3:n+4),fl(2,0:n+1),f(2,0:n+1)
double precision sh(2,n),duml(2,0:n),dumr(2,0:n)
double precision corr(0:n+1),corl(0:n+1),cor(0:n+1)
double precision uap,cap,lap1,lap2,a1,a2,rap1a
double precision rap1b,rap2a,rap2b,b1,b2
double precision w1a,w1b,w2a,w2b,sum1,sum2,g
double precision zm,up,um,u,tmp,c,cm,cp,hap,l1,l2
double precision f1abs,f2abs,tmpl,tmpr,dh,dhu
double precision uas11,uas12,uas13,uas21,uas22,uas23
double precision uas21p,uas22p,uas23p
double precision CCL,CCR,DB,threshold,SL,SR,lap11,lap22
double precision dt,dx,fh(0:n),nutot(n),energyf(n)
double precision slope11, slope12, slope21, slope22
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
        uas13=bed(CL)
        uas21=pv(1,CR)
        uas22=pv(2,CR)
        uas23=bed(CR)
        tmporcl=0
        tmporcr=0
        slope = 3

        if(order.eq.2) then

                !-------------------------slope muscl----------------
                if(slope.eq.1) then !centered slope
                        slope11 = (pv(1,CL+1)-pv(1,CL-1))/(2*dx)
                        slope12 = (pv(2,CL+1)-pv(2,CL-1))/(2*dx)
                        slope21 = (pv(1,CR+1)-pv(1,CR-1))/(2*dx)
                        slope22 = (pv(2,CR+1)-pv(2,CR-1))/(2*dx)
                endif

                if(slope.eq.2) then !Upwind
                        slope11 = (pv(1,CL)-pv(1,CL-1))/(dx)
                        slope12 = (pv(2,CL)-pv(2,CL-1))/(dx)
                        slope21 = (pv(1,CR)-pv(1,CR-1))/(dx)
                        slope22 = (pv(2,CR)-pv(2,CR-1))/(dx)
                endif

                if(slope.eq.3) then !downwind
                        slope11 = (pv(1,CL+1)-pv(1,CL))/(dx)
                        slope12 = (pv(2,CL+1)-pv(2,CL))/(dx)
                        slope21 = (pv(1,CR+1)-pv(1,CR))/(dx)
                        slope22 = (pv(2,CR+1)-pv(2,CR))/(dx)
                endif

                !--------------------------2 order-------------------

                uas11 = pv(1,CL) + slope11 *dx/2
                uas12 = pv(2,CL) + slope12 *dx/2

                uas21 = pv(1,CR) + slope21 * (-dx/2)
                uas22 = pv(2,CR) + slope22 * (-dx/2)      
        
        endif


        if(pv(1,CL).LE.threshold) then 
                uas11=0.0d0!pv(1,CL)   
                uas12=0.0d0!pv(2,CL) ! DRY LEFT
                uas13=bed(CL)
        endif
        if(pv(1,CR).le.threshold) then 
                uas21=0.0d0!pv(1,CR)   
                uas22=0.0d0!pv(2,CR) ! DRY LEFT
                uas23=bed(CR)
        endif

        if(qp(1,CL).le.0.0d0 .and.qp(1,cr).gt.0.0d0) then
                uas11=0.0d0
                uas12=0.0d0
                uas13=bed(CL)
                uas21=pv(1,CR)
                uas22=pv(2,CR)
                uas23=bed(CR)
        endif
        if(qp(1,CR).eq.0.0d0 .and. qp(1,cl).gt.0.0d0) then
                uas21=0.0d0
                uas22=0.0d0
                uas23=bed(CR)
                uas11=pv(1,CL)
                uas12=pv(2,CL)
                uas13=bed(CL)
        endif
        if(iwd(cl).eq.1.and.iwd(cr).eq.1.and.pv(1,cl-1).eq.0.0d0) then
                uas11=pv(1,CL)
                uas12=pv(2,CL)
                uas13=bed(CL)
                tmporcl=1 
        endif
        if(iwd(cl).eq.1.and.iwd(cr).eq.1.and.pv(1,cr+1).eq.0.0d0) then
                uas21=pv(1,CR)
                uas22=pv(2,CR)
                uas23=bed(CR)
                tmporcr=1 
        endif


        if(iwd(cl).eq.1 .and. dabs(uas12).le.1.d-12 .and. tmporcl.eq.0) then 
                uas11=pv(1,cl)-(uas13-bed(cl)) 
        endif
        if(iwd(cr).eq.1 .and. dabs(uas22).le.1.d-12 .and. tmporcr.eq.0) then 
                uas21=pv(1,cr)-(uas23-bed(cr)) 
        endif
!-------------------WET/DRY TREATMENT FOR THE FLUIDS IN MOTION------------
        if(uas11.gt.0.0d0 .and.uas21.eq.0.0d0 .and.uas11.lt.(uas23-uas13)) then
                uas12=0.0d0
                uas22=0.0d0
        endif
        if(uas21.gt.0.0d0 .and. uas11.eq.0.0d0 .and.uas21.lt.(uas13-uas23)) then
                uas12=0.0d0
                uas22=0.0d0
        endif
!--------------------------------------------------------------------------
            
	if(dsqrt(uas11) + dsqrt(uas21) .eq. 0.0d0) then
		uap=0.0d0
	else
		uap = (dsqrt(uas11)*uas12 + dsqrt(uas21)*uas22)/ (dsqrt(uas11) + dsqrt(uas21)) !Computation for F_{I+1/2}
  	endif

  	hap = 0.5d0*(uas11 + uas21)
  	cap = dsqrt(g*hap)
  	lap1 = uap - cap
  	lap2 = uap + cap
  	l1 = uas12 - dsqrt(g*uas11)
  	l2 = uas22 + dsqrt(g*uas21)
  	CCL=dsqrt(g*uas11)  
  	CCR=dsqrt(g*uas21)  
  	SR=dmax1(dmax1(uas22+CCR,lap2),0.0d0)
  	SL=dmin1(dmin1(uas12-CCL,lap1),0.0d0)

  	if(cap.eq.0.0d0) then
   		a1=0.0d0
   		a2=0.0d0
  	else
   		a1 =((lap2*(uas21-uas11)) - (uas22*uas21-uas12*uas11))/(2.0d0*cap)! LeVeque
		a2 =((-lap1*(uas21-uas11)) + (uas22*uas21-uas12*uas11))/(2.0d0*cap)
	endif

	f(1,CL)= uas11*uas12
  	f(2,CL)= (uas11*uas12**2) + 0.5d0*g*(uas11**2)
  	f(1,CR)= uas21*uas22
  	f(2,CR)= (uas21*uas22**2) + 0.5d0*g*(uas21**2)

  	rap1a= 1.0d0
  	rap1b= lap1
  	rap2a= 1.0d0
  	rap2b= lap2
!-------Entropy correction-------	
	if(entropy.eq.3) call hartenl(lap1,lap2,a1,a2,w1a,w1b,w2a,w2b,rap1a,rap1b,rap2a,rap2b,uas22,CCR,uas12,CCR,uap,cap,f1abs,f2abs,i)  
	if(entropy.eq.4)call noentl(lap1,lap2,a1,a2,w1a,w1b,w2a,w2b,rap1a,rap1b,rap2a,rap2b,f1abs,f2abs,i)

	sum1 = w1a + w2a
  	sum2 = w1b + w2b
  	fl(1,i) = (0.50d0*(f(1,CL)+f(1,CR))) - (0.50d0*sum1)
  	fl(2,i) = (0.50d0*(f(2,CL)+f(2,CR))) - (0.50d0*sum2) ! Fluxes

!----------------------------bed fluxes----------------------------------------
! computes the bed fluxes + and - in each face
!
!
!--------------------redefinition of topography--------------------------------
	DB=uas23-uas13

        if(uas11.lt.DB  .and. uas21 .lt. threshold) then
        	b1 =  0.5d0*cap*uas11
          	b2 = -0.5d0*cap*uas11
        elseif(uas21.lt.-DB .and. uas11 .lt. threshold) then
          	b1 = -0.5d0*cap*uas21
          	b2 =  0.5d0*cap*uas21
        else
          	b1 =  0.5d0*cap*DB
          	b2 = -0.5d0*cap*DB
	endif
!------------------------------------------------------------------------------

        duml(1,i) = 0.5d0*(b1*rap1a*(1.0d0-dsign(1.0d0,lap1))+b2*rap2a*(1.0d0-dsign(1.0d0,lap2)))
        duml(2,i) = 0.5d0*(b1*rap1b*(1.0d0-dsign(1.0d0,lap1)) +b2*rap2b*(1.0d0-dsign(1.0d0,lap2)))!R^{-}_{i+1/2}

        dumr(1,i) = 0.5d0*(b1*rap1a*(1.0d0+dsign(1.0d0,lap1))+b2*rap2a*(1.0d0+dsign(1.0d0,lap2)))  !R^{+}_{i+1/2}
        dumr(2,i) = 0.5d0*(b1*rap1b*(1.0d0+dsign(1.0d0,lap1))+b2*rap2b*(1.0d0+dsign(1.0d0,lap2)))
      
ENDDO

DO I=1,N
	SH(1,I) = duml(1,i) +dumr(1,i-1)   !contribution to each cell from bed fluxes
      	SH(2,I) = duml(2,i) +dumr(2,i-1)  
ENDDO

return
end

