module param
	
	integer,parameter :: row=8,crow=8
	integer,parameter :: sca=row*crow
	integer,parameter :: deg=3
end module
include "networklibstandard.f90"
program main
	use param
	use biolibs_fixedscale
	implicit none
	integer :: thenet(sca,deg)
	integer :: i,j,k,mm
	real :: x0,y0,x1,y1,nl,ny,zz
	thenet=0
	do i=0,row-1
		do j=2,crow
			call pedge(thenet,sca,deg,i*crow+j,i*crow+j-1)
		end do
		!period_row:
		call pedge(thenet,sca,deg,i*crow+1,i*crow+crow)
	end do
	do i=1,row
		do j=1,crow
			mm=mod(i*crow+j,sca)

				if(mod(i,2)==1 .and. mod(j,2)==1) then
					call pedge(thenet,sca,deg,mm,(i-1)*crow+j)
				end if
				if(mod(i,2)==0 .and. mod(j,2)==0) then
					call pedge(thenet,sca,deg,mm,(i-1)*crow+j)
				end if

			
		end do
	end do
	


go to 101
	if(mod(row,2)==1) then
		call deledge(thenet,sca,deg,crow,crow-1)
		call deledge(thenet,sca,deg,(row-1)*crow+1,crow*(row-1)+2)
	elseif(mod(row,2)==0) then
		call deledge(thenet,sca,deg,crow,crow-1)
		call deledge(thenet,sca,deg,row*crow,row*crow-1)
	end if

	!call disorderzbound(thenet,row,crow,2,0.1)
101	open(unit=11,file='lattice')
	open(unit=13,file='latview1.csv')
	!write the crystal lattice
	do i=1,sca
		do j=1,deg
			write(11,*) thenet(i,j)!(thenet(i,j),',',j=1,deg)
		end do
	end do


	do i=1,sca
			ny=mod(i,crow)
			if(ny==0) ny=crow
			x0=1.732*ny
			y0=2.0*(i/crow)+(-1)**(i/crow+ny+1)*0.5
		do j=1,deg
		if(thenet(i,j)/=0) then
			nl=mod(thenet(i,j),crow)
			zz=thenet(i,j)/crow
			if(nl==0) then
				nl=crow
				!zz=zz-1
			end  if
			x1=1.732*nl
			y1=2.0*(zz)+(-1.0)**(zz+nl+1)*0.5
			write(13,*) i,',',thenet(i,j),',',x0,',',y0,',',x1,',',y1
		end if
		end do
		
	end do


end program main

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
