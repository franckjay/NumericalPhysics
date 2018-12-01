cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Author: Jay Franck  Email: franckjay@gmail.com
cc Course: Physics 580 Instructor: Chad Kishimoto
cc
cc Program Name: Project 6 - TDSE
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc PROGRAM NOTES:
cc __________________________________________________________________
cc
cc This program solves the Time Dependent Schrodinger Equation using
cc the Crank-Nicolson method for a particle in a box. 
cc
cc SUBROUTINES:
cc __________________________________________________________________
cc
cc matinv.f: Inverts the 2Nx2N matrices in the program 
cc
cc franck_proj6_subs.f: Contains two subroutines. 'bigmat' produces
cc two 2Nx2N matrices using the hamiltonians supplied in the main prog
cc and an identity matrix. 'crank' is the crank-nicolson routine which
cc finds Psi(t+1) from Psi(t) and the big matrices.
cc
cc VARIABLES
cc __________________________________________________________________
cc INPUT.............................................................
cc L: Length of one side of the box
cc N: Number of points for the solution
cc dt = size of the time step for the Crank-Nicolson
cc nsteps: number of time steps to go through
cc width: the width of the original gaussian
cc k = constant for the QHO. Should be small (0.01-0.02)
cc x0: Original center of the Gaussian width
cc
cc OUTPUT............................................................
cc rho: Probability density at each step
cc Sigma: Width of the gaussian at each iteration
cc expect: Expectation value (<x>)at a value of x
cc norm: Normalization constant
cc sum______: The summed value over the number of points in the soln
cc QHOexpect: Classical QHO expectation value at a point
cc anSigma: Analytical value of the width of the Gaussian
cc
cc INTERMEDIATE......................................................
cc H: Hamiltonian of the system
cc KE: Kinetic Energy of the system
cc V: Potential energy of the system
cc dx: Step size of x
cc x: Values of x based on dx over -L to L
cc Big: LHS Big Matrix thats constructed with subroutine bigmat
cc invBig: Inverse of the Big matrix
cc rightBig: RHS Big matrix
cc Psi: Wavefunction vector. Real Part 1 to N, Imag from N+1 to 2N
cc Pi: A piece of a circular Pie
cc dum: dummy vector passed to matinv.f
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	program TDSE

	implicit none
	integer maxsize
	real hbar,m,c
	parameter (maxsize=1000,hbar=1.0,m=1.0,c=1.0)
	

c.......INPUTS........................................................
	real L,dt,x0,width,k
	integer N,nsteps


c.......OUTPUT........................................................

	real rho(maxsize),sigma(maxsize),expect(maxsize),norm(maxsize)
	real sumRho,sumSigma,sumExpect,sumNorm,t,anSigma
	real QHOexpect

c.......INTERMEDIATE..................................................
	real H(maxsize,maxsize),KE(maxsize,maxsize),V(maxsize,maxsize)
	real dx,x(maxsize),Big(maxsize,maxsize),invB(maxsize,maxsize)
	real psi(maxsize),PI,rightBig(maxsize,maxsize)
	integer i,j,dum(maxsize)

	

cc......USER.INPUT.SECTION............................................

	print *,"This is computational project 6"
	print *,"Written by Jay Franck"
	print *,"************************************"
	print *," "


	

c	Input for L
	print *,"How big do you want the box to be?"
	print *,"Enter 'L' (the box will go from -L to L): "
	read *,L
c	Error trap for L
	if (L .LE. 0.0) then
		print *,"Box sides must be positive!"
		GO TO 29
	endif

c	Input of N
	print *,"Enter the number of points for the solution: "
	read *,N

c	Error trap of N
	if (N .LE. 0) then
		print *,"Lattice Points must be positive!"
		GO TO 29
	endif

	if (MOD(N,1) .NE. 0) then
		print *,"Number of lattice points must be an integer!"
		GO TO 29
	endif

c	Input of dt
	print *,"Enter the size of the time step 'dt': "
	read *,dt

c	Error trap of dt
	if (dt .LE. 0) then
		print *,"Time step must be positive!"
		GO TO 29
	endif

c	Input of nsteps
	print *,"Enter the number of time steps: "
	read *,nsteps

c	Error trap of nsteps
	if (nsteps .LE. 0) then
		print *,"Number of time steps must be positive!"
		GO TO 29
	endif

	if (MOD(nsteps,1) .NE. 0) then
		print *,"Number of time steps must be an integer!"
		GO TO 29
	endif

	print *,"The center of the Gaussian is at: "
	read *,x0

	if (x0 .LE. -L .OR. x0 .GE. L) then
		print *,"Center must be between -L and L!"
		GO TO 29
	endif

	print *,"Enter the width of the Gaussian: "
	read *,width
	if (width .LE. 0) then
		print *,"Width must be positive!"
		GO TO 29
	endif

	print *,"Enter a small value for 'k': "
	read *,k
	if (k .LE. 0) then
		print *,"k must be positive!"
		GO TO 29
	endif	


c	Some good preset parameters for testing

c	L= 5.0
c	N=100
c	dt=0.1
c	nsteps=50
c	x0=0.0
c	width=0.5
c	k= 0.01

cc......MAIN.PROGRAM.SECTION...........................................

	PI = 2. * dasin(1.d0)

c	Calculates the step size
	dx = 2.0*L/(N-1)

c	Initialize variables

	sumRho=0.0
	sumNorm=0.0
	sumExpect=0.0
	sumSigma=0.0

	rho = 0.0
	sigma = 0.0 
	expect=0.0
	norm=0.0
	t=0.0

	do i = 1,maxsize
        	do j = 1, maxsize			
	     		KE(i,j) = 0.
	     		V(i,j)=0.
			Big(i,j)=0.0

		enddo
		x(i)=0.0
		psi(i)=0.0
		rho(i)=0.0
		sigma(i)=0.0
		expect(i)=0.0
		norm(i)=0.0

	enddo


c	Fills x(i)/initial psi list
	do i=1,N
		x(i)= -L+(i-1)*dx

		psi(i)=((2.0*PI*width)**0.25)*	
     +		 exp(-((x(i)-x0)**2)/(4.0*width))

	enddo

c	Boundary Conditions for box
c	since Wave vanishes at boundaries
	psi(1)=0.0
	psi(N)=0.0
	

c	Produces the tridiagonal matrix for the discretized KE
	do i = 1, N
		do j = 1, N
			if (i.EQ.j) then	
				KE(i,j) = (hbar**2)/m/(dx**2)
			elseif (i.EQ.(j-1) .OR. i.EQ.(j+1)) then
				KE(i,j) = (-0.5*(hbar**2))/m/(dx**2)
			endif
		enddo
	enddo

c	Calculates Potential energy and gets the Hamiltonian
	do i=1,N
		do j=1,N
			if (i .EQ. j) then

c				V for QHO (Advanced Project)
				V(i,j)=0.5*k*(x(i)**2)

c				Basic Project			
c				V(i,j)= 0.0
			endif
			H(i,j)=(dt/2/hbar)*(V(i,j)+KE(i,j))
		enddo
	enddo
c	Calls the subroutine to make the 2Nx2N matrix
	call bigmat(Big,rightBig,H,N)

c	Inverts the LHS 2Nx2N matrix	
	call matinv(Big,invB,2*N,maxsize,dum)


	
c.......Initial values.....................................
	do i=1,N
		rho(i) = psi(i)**2
		norm(i) =  rho(i)*dx
		expect(i) = x(i)*rho(i)*dx
		sigma(i) = (x(i)**2)*rho(i)*dx
	enddo

	do i=1,N
		
		sumNorm= sumNorm +norm(i)
		sumExpect= sumExpect + expect(i)
		sumSigma= sumSigma + sigma(i)
	enddo

c	Calculates the values at t=0.0
	sumExpect=sumExpect/sumNorm
	sumSigma=(sumSigma/sumNorm)-sumExpect**2
	QHOexpect= x0*cos(t*sqrt(k/m))
	anSigma= width

c*******Print to File Section********************************
	
	print*, 'Normalization written to norm.txt'		
	open(unit=1, file='norm.txt')

	write(1,*) '# created by franck_proj6.f'
	write(1,*) '# TDSE'
	write(1,*) '# Time in First column, Norm values in Second'
	write(1,*) '# t                      Norm'
	write(1,*) t,"      ",sumNorm

	print*, 'Expectation value written to expect.txt'		
	open(unit=2, file='expect.txt')

	write(2,*) '# created by franck_proj6.f'
	write(2,*) '# TDSE'
	write(2,*) '# Time in First column, Expectation in Second'
	write(2,*) '# t                      <X>'
	write(2,*) t,"      ",sumExpect


	print*, 'Variance(Width) written to sigma.txt'		
	open(unit=3, file='sigma.txt')

	write(3,*) '# created by franck_proj6.f'
	write(3,*) '# TDSE'
	write(3,*) '# Time in First column, Sigma values in Second'
	write(3,*) '# t                      Width(Sigma)'
	write(3,*) t,"      ",sumSigma


	print*, 'Analytical Variance(Width) written to anSigma.txt'		
	open(unit=4, file='anSigma.txt')

	write(4,*) '# created by franck_proj6.f'
	write(4,*) '# TDSE'
	write(4,*) '# Time in First column, Sigma values in Second'
	write(4,*) '# t                      Width(Sigma)'
	write(4,*) t,"      ",anSigma

	print*, 'Classic Harmonic Oscillator written to QHOexpect.txt'		
	open(unit=6, file='QHOexpect.txt')

	write(6,*) '# created by franck_proj6.f'
	write(6,*) '# Expectation values for the Harmonic Oscillator'
	write(6,*) '# k = ',k,'                 Initial x = ',x0
	write(6,*) '# Time in First column, <x> values in Second'
	write(6,*) '# t                      QHOexpect'
	write(6,*) t,"      ",QHOexpect
	
c*******End of Print to File Section********************************




c......Crank-Nicolson.......................................

	do i=1,nsteps
		t=t+dt

		call crank(rightBig,invB,N,psi)

c		Initializes values
		sumNorm = 0.0
		sumExpect = 0.0
		sumSigma = 0.0


c		Everything after j = N is the imaginary part
c		Everything before is real in the 2N Matrix
		do j=1,N
			rho(j) = psi(j)**2 + psi(j+N)**2
			norm(j) =  rho(j)*dx
			expect(j) = x(j)*rho(j)*dx
			sigma(j) = (x(j)**2)*rho(j)*dx

		enddo

		do j=1,N
			
			sumNorm= sumNorm + norm(j)
			sumExpect= sumExpect + expect(j)
			sumSigma= sumSigma + sigma(j)
		enddo

c		Calculates the values at this time step		
		sumExpect=sumExpect/sumNorm
		sumSigma=(sumSigma/sumNorm)-sumExpect**2
		ansigma=width +((hbar**4)*(t**2)/(4*width))
		QHOexpect= x0*cos(t*sqrt(k/m))
		
c		Writes calculations to files
		write(1,*) t,"      ",sumNorm
		write(2,*) t,"      ",sumExpect
		write(3,*) t,"      ",sumSigma
		write(4,*) t,"      ",anSigma
		write(6,*) t,"      ",QHOexpect


		

	enddo

c*******End of Crank-Nicolson***************************************

	print*, '# Probability Density written to probDensity.txt'		
	open(unit=5, file='probDensity.txt')

	write(5,*) '# created by franck_proj6.f'
	write(5,*) '# TDSE'
	write(5,*) '# X in First column, Rho values in Second'
	write(5,*) '# x                      Probability Density (Rho)'

c	Prints wavefunction at final step
	do i=1,N
		write(5,*) x(i),"      ",rho(i)
	enddo




c.......OUTPUT...................................................




	close(1)
	close(2)
	close(3)
	close(4)
	close(5)

      	return
29     	end
	




cc......SAMPLE OUTPUT.............


c Phys580/franck_proj6> gfortran franck_proj6.f franck_proj6_subs.f matinv.f
c matinv.f:91.72:

c        IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'                       
c                                                                        1
c Warning: Deleted feature: PAUSE statement at (1)
c Phys580/franck_proj6> ./a.out
c This is computational project 6
c Written by Jay Franck
c ************************************
  
c How big do you want the box to be?
c Enter 'L' (the box will go from -L to L): 
c 5
c Enter the number of points for the solution: 
c 100
c Enter the size of the time step 'dt': 
c  .1
c Enter the number of time steps: 
c 50
c The center of the Gaussian is at: 
c 0.0
c Enter the width of the Gaussian: 
c 0.5
c Enter a small value for 'k': 
c 0.01
c Normalization written to norm.txt
c Expectation value written to expect.txt
c Variance(Width) written to sigma.txt
c Analytical Variance(Width) written to anSigma.txt
c Classic Harmonic Oscillator written to QHOexpect.txt





cc...............................








