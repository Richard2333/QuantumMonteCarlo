module param
integer,parameter :: sca=400
integer,parameter :: leng=20
integer,parameter :: deg=3
end module
include "networklibstandard.f90"
program main
use param
use biolibs_fixedscale
	implicit none
	integer :: thenet(sca,deg)
	integer :: i,j,k
	thenet=0
	do i=0,leng-1
		do j=1,leng-1
			call pedge(thenet,sca,deg,i*leng+j,i*leng+j+1)
		end do
	end do
	do i=0,leng-1
		do j=1,leng-1
			if((i+1)*leng+j<=sca) then
				if(mod(i,2)==0 .and. mod(j,2)==1) then
					call pedge(thenet,sca,deg,i*leng+j,(i+1)*leng+j)
				end if
				if(mod(i,2)==1 .and. mod(j,2)==0) then
					call pedge(thenet,sca,deg,i*leng+j,(i+1)*leng+j)
				end if
			end if
		end do
	end do
	if(mod(leng,2)==1) then
		call deledge(thenet,sca,deg,leng,leng-1)
		call deledge(thenet,sca,deg,(leng-1)*leng+1,leng*(leng-1)+2)
	elseif(mod(leng,2)==0) then
		call deledge(thenet,sca,deg,leng,leng-1)
		call deledge(thenet,sca,deg,leng*leng,leng*leng-1)
	end if

	
	open(unit=11,file='latticenet1.csv')
	!write the crystal lattice
	do i=1,sca
		write(11,*) (thenet(i,j),',',j=1,deg)
	end do



end program main
