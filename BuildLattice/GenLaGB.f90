module param
	integer,parameter :: x=15,y=x
	contains

	function f(x)

	end function f



end module
include "networklibstandard.f90"
program main
	use param
	implicit none
	external sgeev
	integer,parameter :: period=4
	integer :: lattice(x*y,3),exec(y),dx(period),xgb,lend,i,j,ny,nl
	real :: x0,y0,x1,y1
	dx(1)=-1
	dx(2)=-1
	dx(3)=0
	dx(4)=0
	xgb=10
	open(unit=12,file='nets.csv')
	open(unit=13,file='lat.csv')
	lattice=0
	lend=7
	call alattice(lattice,x,y,lend)
	exec=0

	call gengb(lattice,x,y,period,dx,xgb,lend,exec,2)

	!allocate(work(lwork))
	do i=1,x*y
		write(12,*) (lattice(i,j),',',j=1,3)
	end do
	do i=1,x*y
			ny=mod(i,y)
			if(ny==0) ny=y
			x0=1.732*ny
			y0=2.0*(i/y)+(-1)**(i/y+ny+1)*0.5
		do j=1,3
		if(lattice(i,j)/=0) then
			nl=mod(lattice(i,j),y)
			if(nl==0) nl=y
			x1=1.732*nl
			y1=2.0*(lattice(i,j)/y)+(-1.0)**(lattice(i,j)/y+nl+1)*0.5
			write(13,*) i,',',lattice(i,j),',',x0,',',y0,',',x1,',',y1
		end if

	end do
	end do
end


subroutine connect(lattice,x,y,a,b)
	integer ::x,y,lattice(x*y,3),a,b,i,j,k,t
	t=0
	do i=1,3
		if(lattice(a,i)==b) t=1
	end do
	do i=1,3
		if(lattice(a,i)==0 .and. t==0) then
			lattice(a,i)=b
			t=1
		end if
	end do
	t=0
	do i=1,3
		if(lattice(b,i)==a) t=1
	end do
	do i=1,3
		if(lattice(b,i)==0 .and. t==0)then
			lattice(b,i)=a
			t=1
		end if
	end do
	002 end subroutine connect

subroutine disconnect(lattice,x,y,a,b)
	integer :: x,y,lattice(x*y,3),a,b,i,j,k
	do i=1,3
		if(lattice(a,i)==b)then
			lattice(a,i)=0
		end if
	end do

	do j=1,3
		if(lattice(b,j)==a)then
			lattice(b,j)=0
			write(*,*) 'i'
		end if
	end do
	end subroutine disconnect

subroutine removeatom(lattice,x,y,a,b,sft)
	implicit none
	integer :: lattice(x*y,3),x,y,a,b,i,j,k,sft
	if(b>1 .and. b<y) then
		do i=1,3
			k=lattice((a-1)*y+b,i)
			call disconnect(lattice,x,y,(a-1)*y+b,k)
		end do
		call connect(lattice,x,y,(a-1)*y+b+1+sft,(a-1)*y+b-1)

	end if

	end subroutine removeatom

subroutine insertatom(lattice,x,y,a,b)!insert a atom between a,b and a,b+1
	implicit none
	integer :: lattice(x*y,3),x,y,a,b,i,j,k
	call disconnect(lattice,x,y,(a-1)*y+b,(a-1)*y+b+1)

	end subroutine insertatom

	!Magical region following,TAKE CAREFUL
subroutine gengb(lattice,x,y,period,dx,xgb,lend,exec,n)!generate a lattice set with GB
	implicit none

		integer :: v1(2),v2(2),x,y,period,dx(period),xgb
		integer :: lattice(x*y,3),i,j,k,exec(y),lend,sft,n
		sft=0
		do k=1,n
		do i=1,period
			call plus1row(lattice,lend,exec,x,y,mod(lend,2))
			if(dx(i)==-1) then
				call removeatom(lattice,x,y,lend+1,xgb,sft)
				exec(xgb)=1
				xgb=xgb-1
				sft=sft+1
				write(*,*) xgb
			end if
			lend=lend+1
		 end do
		end do

		!sft=sft+1
	end subroutine gengb

subroutine plus1row(lattice,lend,exec,x,y,odd)
	implicit none!exec:index of atom that are excluded
	integer :: lattice(x*y,3),x,y,lend,odd,exec(y),i,j,k
	integer,external :: ifconn
	do i=lend*y+1,lend*y+y
		if(exec(i-y*lend)==1) then
			do j=i,lend*y+y-1
				lattice(j,:)=lattice(j+1,:)
			end do
		end if
		if(exec(i-y*lend)==0) then
			if(odd==1 .and. mod(i-y*lend,2)==1) then
				call connect(lattice,x,y,i,i-y)
			end if
			if(odd==0 .and. mod(i-y*lend,2)==0) then
				call connect(lattice,x,y,i,i-y)
			end if
		end if
	end do
	do i=lend*y+1,lend*y+y
		do j=lend*y+1,lend*y+y
		if(ifconn(lattice,x,y,i-y,j-y)==1) then
				call connect(lattice,x,y,i,j)
			end if
		end do
	end do

	end subroutine plus1row

subroutine alattice(lattice,x,y,lend)
!生成lattice（无缺陷）
	implicit none
	integer :: lattice(x*y,3),x,y,lend,i,j,k
	do i=0,lend-1
		do j=1,y-1
			call connect(lattice,x,y,i*y+j,i*y+j+1)
		end do
	end do
	do i=0,lend-2
		do j=1,y-1
			if(mod(i+1,2)==mod(j,2)) then
				call connect(lattice,x,y,i*y+j,i*y+j+y)
			end if
		end do
	end do

	end subroutine alattice


function ifconn(lattice,x,y,a,b)
	implicit none
	integer :: lattice(x*y,3),x,y,a,b,ifconn,i
	ifconn=0
	do i=1,3
		if(lattice(a,i)==b) ifconn=1
	end do
	return
	end function ifconn
