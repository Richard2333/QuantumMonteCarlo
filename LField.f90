module LFMontecarlo


	contains
	function det(aaa,m)
	    implicit none
	    complex :: aaa(m,m)
	    integer :: m,n,info,i,ipiv,lda
	    external :: dgetrf,cgetrf
	    real :: det

	    !aaa=reshape((/1,3,2,4/),(/2,2/))
	    n=m
	    lda=m
	    det=1.0d0

	    call cgetrf( m, n, aaa, lda, ipiv, info )

	    do i=1,m
	        det=det*aaa(i,i)
	    end do
	    return

	end function det

	function delta(a,b)
		implicit none
		integer :: a,b,delta
		delta=0
		if(a==b) delta=1
		return
	end function delta

	subroutine mchub(del,k,ms,phi,nx,ny,ntm,lattice,vol,tre)
		implicit none
		integer :: nx,ny,ntm,lattice(nx*ny,3)
		integer :: i,j,a,nn,x,y,t1,t2
		real :: phi(ntm*nx*ny),del,k,ms,vol(nx*ny,nx*ny),tre
		complex :: ep,m(nx*ny*ntm,nx*ny*ntm),exp1,m1(nx*ny*ntm,nx*ny*ntm)
		nn=nx*ny
		do i=1,nx*ny*ntm
			do j=1,nx*ny*ntm
				t1=i/nn
				t2=j/nn
				x=mod(i,nn)
				y=mod(j,nn)
				ep=(0,-1)*del*phi(i)
				m(i,j)=delta(x,y)*(delta(t1,t2)-exp(ep)*delta(t1-1,t2))+ms*del*delta(x,y)*delta(t1-1,t2)
				do a=1,3
					m(i,j)=m(i,j)-(k*del*delta(y,lattice(x,a))*delta(t1-1,t2))
				end do
			end do
		end do
		m1=matmul(m,conjg(transpose(m)))
		exp1=0
		do i=1,ntm
			do j=1,nx*ny
				do a=1,nx*ny
					exp1=exp1+phi((i-1)*nn+j)*(1.0/vol(j,a))*phi((i-1)*nn+a)
				end do 
			end  do 
		end do 
		exp1=-0.5*del*exp1
		tre=det(m1,nx*ny*ntm)*exp(exp1)


	
	end subroutine mchub
	!subroutine containing the monte carlo,from a aux.field phi and 
	!standard hubbard hamiltonian and others/see Phys.Rev.B 2014 89,19,195426

	subroutine mcarb


	end subroutine mcarb



end module LFMontecarlo
program main

end program main
