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
		!tre=det(m1,nx*ny*ntm)*exp(exp1)

	end subroutine detmchubb_coh
	!subroutine containing the monte carlo,from a aux.field phi and 
	!standard hubbard hamiltonian and others/see Phys.Rev.B 2014 89,19,195426

	function det(mat,m)!mat会覆盖
		implicit none
		real :: mat(m,m)
		integer :: m,n,info,i,ipiv(m),lda
		external :: sgetrf,cgetrf
		real :: det

	    !mat=reshape((/1,3,2,4/),(/2,2/))
	    n=m
	    lda=m
	    det=1.0d0
	    !write(*,*) m,n
	    call sgetrf( m, n, mat, lda, ipiv, info )
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


		function sInverse(n, a)  result(ra)

		integer::n,lda,ipiv(n),info,lwork

		real ::a(n,n),ra(n,n),work(n)

	  	ra=a

  		lwork=n

  		lda=n

  		call sgetrf(n, n, ra, lda, ipiv, info)

  		if(info/=0) write(0,*) 'Error occured in cgetrf!',info

  		call sgetri(n, ra, lda, ipiv, work, lwork, info)

  		if(info/=0) write(0,*) 'Error occured in cgetri!',info

	endfunction

	function mexp1(mat,dt,n) result(mexp)
		implicit none 
		real :: mexp(n,n),mat(n,n),maz(n,n),dt
		integer :: n,i
		maz=mat*dt
		mexp=0
		!return
		do i=1,n-1
			mexp(i,i)=1.0
		end do
		do i=1,7
			mexp=mexp-maz
			maz=-matmul(maz,mat)
			maz=maz*dt
			!print *,mat
		end do
		return
	end function mexp1
	

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
		!tre=det(m1,nx*ny*ntm)*exp(exp1)


	end subroutine hybmc

!与相干态et萌卡需要的矩阵dim一样
	subroutine insert_parti_mat(Mat,mt,n,t,cl,ro)!分块矩阵加入元素
		implicit none
		integer :: n,t,cl,ro,i,j,s0,s1
		real :: Mat(n*t,n*t),mt(n,n)
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
		real :: Mat(n*nt,n*nt),ham(n,n),k(n,n),v(n,n),b(n,n)
		real :: s(n,nt),dt,sgm     ,thop,lamda
		lamda=2.7*dt/2
		thop=2.7
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

			b=matmul(mexp1(k,dt,n),mexp1(v,dt,n))
			cl=t
			!b=-b
			if(t/=nt) then
				ro=t+1
				b=-b
			else 
				ro=1
				b=0
			end if
			call insert_parti_mat(Mat,b,n,nt,cl,ro)
			!一个分块的组装
		end do


	end subroutine assemb_mat

	subroutine detmcstep(s,n,nt,dt,lattice,det0,func)
	!func 为需要sample的关联函数
	implicit none
		integer :: n,nt,lattice(n,mxdg),i,j,kk,cl,ro,num
		real :: Mat(n*nt,n*nt),func(n*nt,n*nt)
		real :: s(n,nt),sgm,rm,det1,det0,dt,det2
		open(unit=19,file='ll.csv')
		call radm(rm)
		mat=0
		num=int(rm*n*nt)
		s(mod(num,n)+1,num/n)=-s(mod(num,n)+1,num/n)
		write(*,*) '---------------',num
		sgm=1
		Mat=0
		call assemb_mat(Mat,s,n,nt,dt,lattice,sgm)
		sgm=-1
		det1=det(Mat,n*nt)
		write(*,*) "==",det1
		Mat=0
		call assemb_mat(Mat,s,n,nt,dt,lattice,sgm)
		det2=det(Mat,n*nt)
		det1=det1*det2
		
		write(*,*) 'ss',det1,det2
		!if(det1<0) stop"FUCK"
			call radm(rm)

			if(rm<(min(det1/det0,1.0))) then
				det0=det1
			else
				s(mod(num,n)+1,num/n)=-s(mod(num,n)+1,num/n)
				
			end if
		
			
		
	
	end subroutine detmcstep

end module LFMontecarlo
program main
	use LFMontecarlo
	implicit none
	integer ,parameter :: nx=8,ny=8,n=nx*ny,ntm=19,ndot=1000
	integer :: lattice(n,3),i,j,m
	real :: phi(nx*ny,ntm),del,k,ms,vol(nx*ny,nx*ny),bn,dett(ndot)
	real :: tre,grn(n*ntm,n*ntm),gs(n,n)
	call sr1and()
	del=0.03
	k=-2.7
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
	do i=1,ndot
		call detmcstep(phi,n,ntm,del,lattice,dett(i),grn)
		dett(i+1)=dett(i)
		!write(*,*) dett(i)
		
	end do
	grn=grn/ndot
	gs=0
	do i=1,n*ntm
		m=mod(i,n)
		gs(m,m)=gs(m,m)+grn(i,i)
	end do
	do i=1,n
		write(13,*) i/nx,',',mod(i,nx),',',real(gs(i,i))
	end do
	
end program main
