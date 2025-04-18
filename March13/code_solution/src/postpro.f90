subroutine popr(qp,be,dxi,g,t,it,por)
implicit none
integer n,problem,hr,entropy,limiter,i,it,order,bb(n)
double precision qp(4,-3:n+4),be(-3:n+4),plt,cfl,ls
double precision bp(-3:n+4),an(n),n1,n2,n3,pv(4,-3:n+4),fh(0:n)
double precision dt,dx,a,b,t,tz,x,g,dxi(0:n+1),PHI(n),W(n),R(n)
double precision nutot(n),por(-3:n+4)
character*30 fileo,fout
common /sc/ cfl,tz,plt,ls,order,entropy,limiter
common /dn/ dt,dx,n
common /bd/ a,b

call primitives(qp,pv,por,be)

if(t.ne.tz ) then 
  	fileo='../output/depth'
  	FOUT=TRIM(TRIM(FILEO)//CHAR(48+INT(IT/10))//CHAR(48+MOD(IT,10))//'.m')

 	OPEN(UNIT=101,FILE=FOUT)
 
 	write(101,*) 'eta=['
 	do i=1,n
    		write(101,*) qp(1,i)/por(i)
 	enddo
 	write(101,*) '];'

 	write(101,*) 'primitives=['
 	x=a
 	do i=1,n
    		write(101,*) pv(1,i),pv(2,i)
 		x=x+dx
 	enddo
 	write(101,*) '];'


 	write(101,*) 't=['
 	write(101,*) t
 	write(101,*) '];'
 
 	close(101)

else
 

        write(14,*) 'x=['
        write(11,*) 'eta=['
        write(13,*) 'b=['
        write(15,*) 'por=['

        call primitives(qp,pv,por,be)
        x=a
        do i=1,n
    	        write(14,*) dxi(i)
                write(11,*) qp(1,i)/qp(3,i)
                write(13,*) be(i)
                write(15,*) por(i),qp(3,i)
                write(32,*) x/1.0d0,(qp(1,i)-1.0d0)/1.0d0
                x=x+dx
        enddo
        write(11,*) '];'
        write(11,*) 'primitives=['
        do i=1,n
             write(11,*) pv(1,i),pv(2,i)
        enddo
        write(11,*) '];'

        write(14,*) '];'
        write(13,*) '];'
        write(15,*) '];'

        write(11,*) 't=['
        write(11,*) t
        write(11,*) '];'

endif

return
end
