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

!####################################Implement the modified DOT solver
return
end

