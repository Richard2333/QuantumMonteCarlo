module LFMontecarlo
	integer,parameter :: mxdg=3

contains
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
	
	function mexp(mat,m,dt)
		implicit none
		integer :: m,i
		complex :: mat(m,m),mexp(m,m)
		real :: dt		
		do i=1,m
			mexp(i,i)=1
		end do
		mexp=mexp-mat*dt
	!	mexp=matmul(mexp,mexp)
	!	mexp=matmul(mexp,mexp)
		return
	end function mexp
	
		
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
	    !write(*,*) m,n
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

	function cInverse(n, a)  result(ra)

		integer::n,lda,ipiv(n),info,lwork

		complex ::a(n,n),ra(n,n),work(n)

	  	ra=a

  		lwork=n

  		lda=n

  		call cgetrf(n, n, ra, lda, ipiv, info)

  		if(info/=0) write(0,*) 'Error occured in cgetrf!',info

  		call cgetri(n, ra, lda, ipiv, work, lwork, info)

  		if(info/=0) write(0,*) 'Error occured in cgetri!',info

	endfunction

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

	subroutine assemb_mat(Mat,s,n,nt,dt,lattice,sgm)!组装矩阵M
		implicit none
		integer :: n,nt,t,lattice(n,mxdg),i,j,kk,cl,ro!,mxdg
		complex :: Mat(n*nt,n*nt),ham(n,n),k(n,n),v(n,n),b(n,n)
		real :: s(n,nt),dt,sgm     ,thop,lamda
		thop=2.7
		lamda=2*thop*dt/2
		
		mat=0
		do i=1,n*nt
			mat(i,i)=1.0
		end do
		do t=1,nt
			k=0
			do i=1,n
				do j=1,mxdg
					if(lattice(i,j)/=0) then
						k(i,lattice(i,j))=-thop
						k(lattice(i,j),i)=-thop
					end if
				end do
			end do 
			v=0
			b=0
			do i=1,n
				v(i,i)=(1/dt)*lamda*sgm*s(i,t)
				b(i,i)=(1,0)
			end do 
			
			b=b-dt*k-dt*v+dt*dt*matmul(k,v)/2
			cl=t
			if(t/=nt) then
				ro=t+1
			else 
				ro=1
			end if
			call insert_parti_mat(Mat,b,n,nt,cl,ro)
			!一个分块的组装
		end do


	end subroutine assemb_mat

	subroutine detmcflip(mat,s,n,nt,dt,lattice,det0,wht,num)
	!func 为需要sample的关联函数
	implicit none
		integer :: n,nt,lattice(n,mxdg),i,j,kk,cl,ro,num
		complex :: Mat(n*nt,n*nt),func(n*nt,n*nt)
		real :: s(n,nt),sgm,rm,det1,det0,dt,wht
		open(unit=19,file='ll.csv')

		s(mod(num,n)+1,num/n)=-s(mod(num,n)+1,num/n)
		sgm=1
		call assemb_mat(Mat,s,n,nt,dt,lattice,sgm)
		det1=det(mat,n*nt)
		func=func+cInverse(n*nt,Mat)
		sgm=-1
		
		call assemb_mat(Mat,s,n,nt,dt,lattice,sgm)
		det1=det1*det(mat,n*nt)
		det1=det1/1000000

		wht=det1/abs(det1)
		!det1=abs(det1)
		write(*,*) 'ss',det1
			call radm(rm)

			if(rm<(min(det1/det0,1.0))) then
				det0=det1
			else
				s(mod(num,n)+1,num/n)=-s(mod(num,n)+1,num/n)
				
			end if
		
			
		
	
	end subroutine detmcflip
	
	
	
end module LFMontecarlo
!运行速度超级丧病慢……
program main
	use LFMontecarlo
	implicit none
	integer ,parameter :: nx=8,ny=8,n=nx*ny,ntm=18,ndot=400 
	integer :: lattice(n,3),i,j,m,x0,y0,pp
	real :: phi(nx*ny,ntm),del,k,ms,vol(nx*ny,nx*ny),bn,dett(ndot),wht,wt
	complex :: grn(n*ntm,n*ntm),mat(n*ntm,n*ntm),gs(n,n)
	call sr1and()
	k=-2.7
	del=0.029
	
	ms=0.0
	open(unit=12, file='lattice')
	open(unit=13,file='gnn.csv')
	do i=1,n
		do j=1,3
			read(12,*) lattice(i,j)
		end do 
	end do 
	do i=1,ntm
		do j=1,n
			call radm(bn)
				phi(j,i)=abs(bn-0.5)/(bn-0.5)
		!随机生成初始phi
				!write(*,*) phi(j,i)
		end do
	end do

	! stop "s"
	dett(1)=0
	grn=0
	wt=0
	do i=1,ndot
		do j=1,n
			call detmcflip(mat,phi,n,ntm,del,lattice,dett(i),wht,n)
		end do
		dett(i+1)=dett(i)
		wt=wt+wht
		write(*,*) i,dett(i),wt
		
		grn=grn+cInverse(n*ntm,mat)*dett(i)
	end do
	grn=grn/(sum(dett)*wt/ndot)
	
	gs=0
	do i=1,n
		!m=mod(i,n)
		gs(i,i)=grn(i,i)
	end do
	do i=1,n
		pp=mod(i,ny)
		if(pp==0) pp=ny
		x0=1.732*pp
		m=i/ny
		if(pp==ny) m=m-1 
		y0=2.0*(m)+(-1)**(m+pp+1)*0.5
		write(13,*) x0,',',y0,',',real(grn(i,i))
	end do
	
end program main
