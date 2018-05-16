module LFMontecarlo

contains
	function det(mat,m)
		implicit none
		complex :: mat(m,m)
		integer :: m,n,info,i,ipiv(m),lda
		external :: dgetrf,cgetrf
		real :: det

	    !mat=reshape((/1,3,2,4/),(/2,2/))
	    n=m
	    lda=m
	    det=1.0d0
	    write(*,*) m,n
	    call cgetrf( m, n, mat, lda, ipiv, info )
	    write(*,*) info
	    do i=1,m
	    	det=det*mat(i,i)
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

	subroutine detmchubb_coh(del,k,ms,phi,nx,ny,ntm,lattice,vol,tre)
		!time slice/hopp energy/
		implicit none
		integer :: nx,ny,ntm,lattice(nx*ny,3)
		integer :: i,j,a,nn,x,y,t1,t2
		real :: phi(ntm*nx*ny),del,k,ms,vol(nx*ny,nx*ny)
		complex :: ep,m(nx*ny*ntm,nx*ny*ntm),exp1,m1(nx*ny*ntm,nx*ny*ntm),tre
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
		write(*,*) exp1
		tre=det(m1,nx*ny*ntm)*exp(exp1)

	end subroutine detmchubb_coh
	!subroutine containing the monte carlo,from a aux.field phi and 
	!standard hubbard hamiltonian and others/see Phys.Rev.B 2014 89,19,195426


	subroutine hybmc(del,k,ms,phi,nx,ny,ntm,lattice,vol,tre)
		!time slice/hopp energy/
		implicit none
		integer :: nx,ny,ntm,lattice(nx*ny,3)
		integer :: i,j,a,nn,x,y,t1,t2
		real :: phi(ntm*nx*ny),del,k,ms,vol(nx*ny,nx*ny)
		complex :: ep,m(nx*ny*ntm,nx*ny*ntm),exp1,m1(nx*ny*ntm,nx*ny*ntm),tre
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
		write(*,*) exp1
		tre=det(m1,nx*ny*ntm)*exp(exp1)


	end subroutine hybmc

!与相干态et萌卡需要的矩阵dim一样
	subroutine insert_parti_mat(Mat,mt,n,t,cl,ro)!分块矩阵加入元素
		implicit none
		integer :: n,t,cl,ro,i,j,s0,s1
		complex :: Mat(n*t,n*t),mt(n,n)
		s0=(ro-1)*n
		s1=(cl-1)*n
		do i=1,n
			do j=1,n 
				Mat(s0+i,s1+j)=mt(i,j)
			end do 
		end do

	end subroutine insert_parti_mat

	subroutine assemb_mat(Mat,s,n,t,dt,lattice,mxdg,sgm)!组装矩阵M
		implicit none
		integer :: n,t,lattice(n,mxdg),i,j,kk,mxdg,cl,ro
		complex :: Mat(n*t,n*t),ham(n,n),k(n,n),v(n,n),b(n,n)
		real :: s(n,t),dt,sgm     ,thop,lamda
	
		
		k=0
		do i=1,n
			do j=1,mxdg
				if(lattice(i,j)/=0) k(i,lattice(i,j))=-thop
			end do
		end do 
		v=0
		b=0
		do i=1,n
			v(i,i)=(1/dt)*lamda*sgm*s(i,t)
			b(i,i)=1
		end do 
		
		b=b-dt*k-dt*v!dt*dt*matmul(k,v)
		call insert_parti_mat(Mat,b,n,t,cl,ro)
		!一个分块的组装


	end subroutine assemb_mat



end module LFMontecarlo
program main
	use LFMontecarlo
	implicit none
	integer ,parameter :: nx=40,ny=8,n=nx*ny,ntm=4
	integer :: lattice(n,3),i,j
	real :: phi(ntm*nx*ny),del,k,ms,vol(nx*ny,nx*ny),bn
	complex :: tre
	call sr1and()
	open(unit=12, file='lattice')
	do i=1,n
		do j=1,3
			read(12,*) lattice(i,j)
		end do 
	end do 
	call radm(bn)
	write(*,*) bn
	vol=10000
	do i=1,n
		!do j=1,n
			!call radm(bn)
			vol(i,i)=9.3
		!end do 
	end do 
	do i=1,ntm*n
		call radm(bn)
		phi(i)=bn
	end do 
	del=0.1
	k=2.7
	ms=0.0
	write(*,*) del,k,ms
	call detmchubb_coh(del,k,ms,phi,nx,ny,ntm,lattice,vol,tre)
	write(*,*) tre 
	

end program main
