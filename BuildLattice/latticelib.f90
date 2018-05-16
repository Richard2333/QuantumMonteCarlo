module param

	integer,parameter :: row=8,crow=40
	integer,parameter :: sca=row*crow
	integer,parameter :: deg=3
end module
include "networklibstandard.f90"

subroutine disorderzbound(thenet,row,crow,zl,ds)
	use biolibs_fixedscale
	implicit none
	integer :: thenet(row*crow,3)
	integer :: row,crow,zl,i,j,k
	real :: ds,s
	call sr1and()
	do i=0,zl-1
		do j=1,crow
			call radm(s)
			if(s<=ds) then
				call remove(thenet,row,crow,i*crow+j)
				write(*,*) i*crow+j
			end if
		end do
	end do

	do i=row-zl,row-1
		do j=1,crow
			call radm(s)
			if(s<=ds) then
				call remove(thenet,row,crow,i*crow+j)
				write(*,*) i*crow+j
			end if
		end do
	end do
end subroutine disorderzbound

subroutine remove(thenet,row,crow,n)
	use biolibs_fixedscale
	implicit none
	integer :: thenet(row*crow,3)
	integer :: row,crow,n1,n2,i,j,maxdeg,n
	maxdeg=3
	n1=n
	do i=1,3
		n2=thenet(n,i)
		if(n2/=0) then

			call deledge(thenet,row*crow,maxdeg,n1,n2)
		end if
	end do
end subroutine remove

subroutine
