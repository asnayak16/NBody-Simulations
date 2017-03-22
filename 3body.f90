program threebody
	real(kind=4), parameter :: PI = 4*atan(1.), G = 4.498279e-8
	integer,parameter :: npar = 1000,dim = 3
	real(kind=4),dimension(2) :: M
	real(kind=4), dimension(npar,dim) :: x,v
	integer :: nt, dsnap, np1,np2
	real :: t,dt
	! Chunk to be divided between threads
	integer,parameter :: chunk = 1000

!	Input Data
	call input_data(x,v,np1,np2,M)
	v(3:np1+2,2)	=	v(3:np1+2,2) + v(1,2)
	v(np1+3:np1+np2+2,2) =	v(np1+3:np1+np2+2,2) + v(2,2) 
	
!   Simulation Parameters
    nt = 5000
    dt = 0.1
    t = 0
	dsnap = 20

!   Time-Marching Loop
	do n = 1,nt
		if (modulo(n-1,dsnap)==0) then
			call output(x,v,t)
		endif
		
		!$OMP PARALLEL SHARED(x,v) PRIVATE(i,a,r)
  		!$OMP DO SCHEDULE(STATIC,CHUNK)
		do i=1,2+np1+np2
			call leapfrog(x,v,dt)
		enddo
		!$OMP ENDDO
		!$OMP END PARALLEL
		
		t = t+dt
	end do	

	call output(x,v,t) ! Output last values

contains
!***********************************************************************
subroutine input_data(x,v,np1,np2,M)
! Initialize Position and Velocity values
	real(kind=4), dimension(npar,dim), intent(out) :: x,v
	real(kind=4), dimension(2),intent(out) :: M
	integer, intent(out) :: np1,np2
	integer :: k 

	open(312,FILE='body.dat',STATUS="OLD",FORM="FORMATTED",ACTION="READ")
	read(312,*) x(1,:),v(1,:),M(1)
	read(312,*) x(2,:),v(2,:),M(2)
	close(312)

	open(312,FILE='A.dat',STATUS="OLD",FORM="FORMATTED",ACTION="READ")
	read(312,*) np1
	do k=1,np1 
		read(312,*) x(k+2,:),v(k+2,:)
	enddo
	close(312)

	open(312,FILE='B.dat',STATUS="OLD",FORM="FORMATTED",ACTION="READ")
		read(312,*) np2
	do k = 1,np2 
		read(312,*) x(k+2+np1,:), v(k+2+np1,:)
	enddo
	close(312)

endsubroutine input_data
!***********************************************************************
subroutine leapfrog (x,v,dt)
! Leapfrog using Position-Verlet algorithm
    real(kind=4), dimension(npar,dim),intent(inout) :: x,v
    real(kind=4), intent(in) :: dt
    real(kind=4),dimension(dim) :: a

	x(i,:) = x(i,:) + 0.5 * dt * v(i,:) 
	call accel(x,a)
	v(i,:) = v(i,:) + dt * G * a
	x(i,:) = x(i,:) + 0.5 * dt * v(i,:)
	
endsubroutine leapfrog
!***********************************************************************
subroutine accel(x,a)
! Calculate Gravitation Acceleration based on position
    real(kind=4), dimension(npar,dim),intent(in) :: x
    real(kind=4), dimension(dim),intent(out) :: a
	real(kind=4) :: r
	integer :: j
	
	! Acceleration  
	a = (/0.,0.,0./)
	do j = 1,2
		if (i/=j) then
			r =	 norm2(x(i,:)-x(j,:)) 		 
			a = a - M(j)*(x(i,:)-x(j,:))/(r**3)
		endif
	enddo
   
endsubroutine accel
!***********************************************************************
subroutine output(x,v,t)
    real(kind=4),dimension(npar,dim),intent(in) :: x,v
    real(kind=4),intent(in) :: t
    integer :: k, f_handle = 314
    character(len=12) :: f_name, num_sfx
	
    write(num_sfx,"(I5)") n
	f_name = 'mc_'//trim(adjustl(num_sfx))//'.dat'
    open(f_handle,FILE=trim(f_name),FORM="formatted")
    do k = 1,np1+np2+2
		write(f_handle,FMT="(6E20.6)") x(k,:),v(k,:)
    enddo
	close(f_handle)

endsubroutine output
!***********************************************************************
endprogram threebody
