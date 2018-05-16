!独立粒子近似，基空间为{|i>}
module mod1

contains

function f(x)

end function f

	function dens(evec,n)
	implicit none
	complex :: evec(n)
	real :: dens(n)
	integer :: n
	
	end function dens


end module

program main
use param
	implicit none
	external sgeev

	
	open(unit=12,file='lattice')
	001 call gengb
	
	
	
	end
	
	subroutine genHa(latce,ham,n,t)
	implicit none
	integer :: latce(n,3),n,i,j,k
	complex :: ham
	real :: t
	do i=1,n
		do j=1,3
			if(latce(i,j)/=0) then
				ham(i,latce(i,j))=t
			end if
		end do
	end do
	
	end subroutine genHa
	
	subroutine ldens(evec,evalue,ef,n)
	use mod1
	implicit none
	integer :: i,j,k,n
	real :: ef,ldos(n),evalue(n)
	complex :: evec(n,n)
	do i=1,n
	 if(evalue(i)<=ef) then
	 	ldos=ldos+dens(evec(i,:),n)
	 end if
	end do
	end subroutine ldens
	
	subroutine connect(lattice,x,y,a,b)
	integer ::x,y,lattice(x*y,3),a,b,i,j,k
	do i=1,3
		if(lattice(a,i)==0) then
			lattice(a,i)=b
		end if
		if(lattice(b,i)==0)then
			lattice(b,i)=1
		end if
	end do
	end subroutine connect
	
	subroutine disconnect(lattice,x,y,a,b)
	integer :: x,y,lattice(x*y,3),a,b,i,j,k
	do i=1,3
		if(lattice(a,i)==b)then
			lattice(a,i)=0
		end if
		if(lattice(b,i)==1)then
			lattice(b,i)=0
		end if
	end do
	end subroutine disconnect
	
	subroutine removeatom(lattice,x,y,a,b)
	integer :: lattice(x*y,3),x,y,a,b
	if(a>1 .and. a<x) then
		call disconnect(lattice,x,y,(a-1)*y+b,(a-1)*y+b-1)
		call disconnect(lattice,x,y,(a-1)*y+b,(a-1)*y+b+1)
		call connect(lattice,x,y,(a-1)*y+b+1,(a-1)*y+b-1)
	end subroutine removeatom
	
	subroutine insertatom
	
	end subroutine insertatom
	
	
	!Magical region following,DONOT TOUCH
	subroutine gengb(lattice1,latticegb,x,y,period,dx,x_gb)!generate a lattice set with GB
	implicit none
		integer :: v1(2),v2(2),x,y,period,dx(period),x_gb
		integer :: lattice1(x*y,3),latticegb(x*y,3)
		
	end subroutine gengb
		
	subroutine trans_green(lattice,J)
		
	end subroutine trans_green

!有问题，不完整/缺少自动生成中间键的部分：
	subroutine latticeshift(lattice,v,x,y)
	implicit none
	integer :: lattice(x*y,3),x,y,v(2),i,j,k,m,n,sf(2)
	sf(1)=v(1)-v(2)
	sf(2)=v(1)+v(2)
	do i=2,x
		if(mod(i-1,sf(2))==0) then
			do j=(i-2)*x+1,(i-2)*x+y
				do k=1,3
					if(lattice(j,k)/x==i-1) then
						m=lattice(j,k)
						n=mod(lattice(j,k),y)
						call disconnect(lattice,x,y,lattice(j,k),j)
						if(n+sf(1)<=y .and. n+sf(1)>=1) then 
							call connect(lattice,x,y,j,m-v(1))
						end if
					end if
				end do
			end do
		end if
	end do
	end subroutine latticeshift

	
