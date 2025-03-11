subroutine hartenl(lap1,lap2,a1,a2,w1a,w1b,w2a,w2b,rap1a,rap1b,rap2a,rap2b,u,c,um,cm,uap,cap,f1abs,f2abs,i)
implicit none
integer n,i,dime
double precision lap1,lap2,w1a,w1b,w2a,w2b
double precision phi,delta,a1,a2,rap1a,rap1b,rap2a,rap2b
double precision uap,cap,u,c,um,cm,up,cp
double precision f1abs,f2abs
character side


delta=dmax1(0.0d0,uap-cap-um+cm,u-c-uap+cap)

if(dabs(lap1).ge.delta) then
  	phi = dabs(lap1)
else
  	phi = (lap1**2 + delta**2)/(2.0d0*delta)
endif
f1abs=phi
w1a = phi*a1*rap1a
w1b = phi*a1*rap1b

delta=dmax1(0.0d0,uap+cap-um-cm,u+c-uap-cap)

if(dabs(lap2).ge.delta) then
 	 phi = dabs(lap2)
else
  	phi = (lap2**2 + delta**2)/(2.0d0*delta)
endif
f2abs=phi
w2a = phi*a2*rap2a
w2b = phi*a2*rap2b

return
end
!------------------------------------------------------------------------
subroutine noentl(lap1,lap2,a1,a2,w1a,w1b,w2a,w2b,rap1a,rap1b,rap2a,rap2b,f1abs,f2abs,i)
implicit none
integer n,i,dime
double precision lap1,lap2,w1a,w1b,w2a,w2b
double precision rap1a,rap1b,rap2a,rap2b
double precision a1,a2,f1abs,f2abs

w1a =dabs(lap1)*a1*rap1a
w1b =dabs(lap1)*a1*rap1b
w2a =dabs(lap2)*a2*rap2a
w2b =dabs(lap2)*a2*rap2b

f1abs = dabs(lap1)
f2abs = dabs(lap2)

return
end
