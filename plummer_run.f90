! Program to simulate Plummer's sphere (Leapfrog integrator)
!! Written by Ashwin Nayak, asnayak@ucsd.edu
!! Compile with -fdefault-real-8 (gfortran) or -real_size=64 (icc)
!! Last Update: 25-feb-2017
! ---------------------------------------------------------
!! OpenMP instructions : 
!! Set env variable OMP_NUM_THREADS=8 (or accordingly)
!! Set env variable OMP_STACKSIZE=2000000 
!! Compile with -fopenmp (gfortran) or -Qopenmp (ifort) switch 
! ---------------------------------------------------------

program plummer
    IMPLICIT NONE
    real(kind=8), parameter :: PI = 4*atan(1.)
	real(kind=8), parameter :: G = 4.498279e-3
    integer, parameter :: dim = 3, npar=10000
    real(kind=8), dimension(npar,dim) :: x,v
	real(kind=8) :: M_tot = 1e11, m 
    integer :: i,j,n, nt, dsnap
    real(kind=8) :: r,t, dt
	real, dimension(dim) :: a  
	
	! Chunk to be divided between threads
	integer,parameter :: chunk = 1000
	
	! Input File
	character(len=12) :: filename = 'pl0.dat'

!   Setup Initial Positions and Velocities and masses
	open(312,FILE=trim(filename),STATUS="OLD",FORM="FORMATTED",ACTION="READ")
	do j = 1,NPAR
		read(312,*) x(j,:)
		read(312,*) v(j,:)
	enddo
	close(312)
	!v = v*0.
    m = M_tot/npar
	
!   Simulation Parameters
    nt = 3000
    dt = 0.005
    t = 0.
	dsnap = 3

!   Time-Marching Loop
    do n = 1,nt
		!$OMP PARALLEL SHARED(x,v) PRIVATE(i,a,r)
  		!$OMP DO SCHEDULE(STATIC,CHUNK)
		do i=1,npar
		! ---------------------------------------------------------
		!				Position Verlet Algorithm 
			x(i,:) = x(i,:) + 0.5 * dt * v(i,:) 
		
		!  Acceleration due to other masses
			a = (/0.,0.,0./)
			do j=1,npar
				if (i/=j) then
					r =	 norm2(x(i,:)-x(j,:)) 		 
					a = a - (x(i,:)-x(j,:))/(r**3)
				endif
			enddo
		
			v(i,:) = v(i,:) + dt * G * m * a
			x(i,:) = x(i,:) + 0.5 * dt * v(i,:)
		! ---------------------------------------------------------
		enddo
		!$OMP ENDDO
		!$OMP END PARALLEL
		t = t+dt
		
		if (modulo(n,dsnap)==0) then
            call output(x,t)
			if(modulo(n,90)==0) then
				print*, 'Iteration: ', n
			endif
        endif
    enddo

    call output(x,t) ! Output last values

contains
!***********************************************************************
subroutine output(x,t)
    real(kind=8),dimension(npar,dim),intent(in) :: x
    real(kind=8),intent(in) :: t
    integer :: k, f_handle = 314
    character(len=12) :: f_name, num_sfx
!    do k = 1,npar	
!		write(*,*) k,t,x(k,:)
!	enddo	
    write(num_sfx,"(I5)") n
	f_name = 'pl'//trim(adjustl(num_sfx))//'.dat'
    open(f_handle,FILE=trim(f_name),FORM="formatted")
    do k = 1,npar
		write(f_handle,FMT=*) x(k,:)
		write(f_handle,FMT=*) v(k,:)
    enddo
	close(f_handle)
endsubroutine output
!***********************************************************************
endprogram plummer
