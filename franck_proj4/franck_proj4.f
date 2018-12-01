cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Author: Jay Franck  Email: franckjay@gmail.com
cc Course: Physics 580 Instructor: Chad Kishimoto
cc
cc Program Name: Project 4 - Schrodinger-Where's the Cat?
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc PROGRAM NOTES:
cc __________________________________________________________________
cc 
cc Calculates the energy levels of orbits for the Woods-Saxon
cc Potential.
cc 
cc SUBROUTINES:
cc __________________________________________________________________
cc 
cc eiglib.f : Finds the eigenvalues and vectors of a matrix and sorts
cc  the values using eigsort from min to max.
cc
cc VARIABLES
cc __________________________________________________________________
cc INPUT.............................................................
cc 
cc N : Number of points for the solution
cc L : The radius of the square well with center at x=0
cc R : Radius of the nucleus in units of fm
cc 
cc OUTPUT............................................................
cc
cc energ(#) : List of the energies of the system, where (1) is the 
cc  ground state, (2) is first excited, etc.
cc evectors(#,#): Gives the e-vectors for each state, where ev(#,1)
cc  represents the ground state, ev(#,2) for the first excited, etc.
cc
cc
cc INTERMEDIATE......................................................
cc 
cc x(#): Represents the position along the bottom of the well, x(1)
cc  is the left edge of the well (position at -L)
cc KE(#,#): Kinetic Energy Tri-Diagonal Matrix
cc V(#,#): Potential Energy Diagonal matrix
cc Ham (#,#) : The Hamiltonian tri-diagonal operator that is sent to
cc  eiglib.f to find the evalues and vectors
cc work(#) : dummy matrix
cc integers i,j,rad: Iteration integers
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	program schrodinger

	implicit none
	integer maxsize
	parameter (maxsize=1000)
c.......INPUTS........................................................
	integer N
	real L,R

c.......OUTPUT........................................................
	real evector(maxsize,maxsize),energ(maxsize)

c.......INTERMEDIATE..................................................
	real PI
	real hbar,Vo,m,a,k,dx
	parameter(hbar = 197.33) 
	parameter(Vo = 50.0)
	parameter(m = 939.0) 
	parameter(a = 0.2) 
	parameter(k = (hbar**2)/m)

	real x(maxsize),Ham(maxsize,maxsize),V(maxsize,maxsize)
	real KE(maxsize,maxsize),work(maxsize)

	integer i,j,rad

 

cc......USER.INPUT.SECTION............................................

	print *,"This is computational project 4"
	print *,"Written by Jay Franck"
	PI = 2. * dasin(1.d0)


c	Input for L
	print *,"How big do you want the box to be?"
	print *,"Enter 'L' (the box will go from -L to L): "
	read *,L
c	Error trap for L
	if (L .LE. 0.0) then
		print *,"Box sides must be positive!"
		GO TO 29
	endif
c	Input of R
	print *,"Enter the nuclear radius (R) in units of fm: "
	read *,R
c	Error trap of R
	if (R .LE. 0.0) then
		print *,"Radius must be positive!"
		GO TO 29
	endif

c	Input of N
	print *,"Enter the number of points for the solution: "
	read *,N
c	Error trap of N
	if (R .LE. 0) then
		print *,"Radius must be positive!"
		GO TO 29
	endif

	if (MOD(N,1) .NE. 0) then
		print *,"Radius must be an integer!"
		GO TO 29
	endif
	

cc......MAIN.PROGRAM.SECTION...........................................



c	Calculates the step size
	dx = 2*L/(N-1)


c	Fills x(i) list
	do i=1,N
		x(i)= -L+(i-1)*dx

	enddo

c	Initializes the variables
	do i = 1,maxsize
        	do j = 1, maxsize
	     		KE(i,j) = 0.
	     		V(i,j)=0.
			energ(i)=0.
			evector(i,j)=0.
			work(i)=0.
		enddo
	enddo

c...............Kinetic Energy

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


c...............Potenital Energy

	open(unit=7, file='ground.dat')
	open(unit=8, file='first.dat')
	open(unit=9, file='second.dat')
	open(unit=10, file='third.dat')
	
	write(7,*) "# Produced by franck_proj4.f"
	write(7,*) "L= ",L,"N: ",N
	write(7,*) "# Radius    Energy"
	write(7,*) "###################"
	write(8,*) "# Produced by franck_proj4.f"
	write(8,*) "L= ",L,"N: ",N
	write(8,*) "# Radius    Energy"
	write(8,*) "###################"
	write(9,*) "# Produced by franck_proj4.f"
	write(9,*) "L= ",L,"N: ",N
	write(9,*) "# Radius    Energy"
	write(9,*) "###################"
	write(10,*) "# Produced by franck_proj4.f"
	write(10,*) "L= ",L,"N: ",N
	write(10,*) "# Radius    Energy"
	write(10,*) "###################"



c	Calculates the potential energy, Hamiltonian, finds the 
c	eigenvalues and evectors, sorts them in order, and writes
c	the values to their respective files for a given nuclear
c	radius

	do rad=2,18
		do i=1,N
			do j=1,N
				if (i .EQ. j) then

c				V for standard QHO
c				V(i,j)=0.5*k*(x(i)**2)
				
				V(i,j)= -Vo/(1+exp((abs(x(i))-R)/a))
				endif
				Ham(i,j)=V(i,j)+KE(i,j)
			enddo
		enddo

		call eig(Ham,N,maxsize,energ,evector,work)
		call eigsrt(energ,evector,N,maxsize)


		write(7,*) R," ",energ(1)
		write(8,*) R," ",energ(2)
		write(9,*) R," ",energ(3)
		write(10,*) R," ",energ(4)
c		Increase the radius
		R=R+0.5
	enddo

	
	print *,"Data written to: 'ground.dat', 'first.dat',"
	print *," 'second.dat', and 'third.dat'"

	
	
c.......OUTPUT.........................................................

c	For Testing Purposes	

c	print *, "Ground State: ", (hbar**2)*(PI**2)/8.0/m/(L**2)
c	print *, "First: ", 4.0*(hbar**2)*(PI**2)/8.0/m/(L**2)

c	print *, "Ground State: ", 0.5*(hbar**2)/m
c	print *, "First: ", 1.5*(hbar**2)/m
	print *,"The first four energy levels for R= ",(R-0.5)
	
	print *,"Ground State Energy (MeV)= ",energ(1)
	print *,"First Excited Energy (MeV)= ",energ(2)	
	print *,"Second Excited Energy (MeV)= ",energ(3)
	print *,"Third Excited Energy (MeV)= ",energ(4)

	close(7)
	close(8)
	close(9)
	close(10)
      	return
29     	end




cc......SAMPLE OUTPUT.............

c Phys580/franck_proj4> gfortran franck_proj4.f eiglib.f
c Phys580/franck_proj4> ./a.out
c This is computational project 4
c Written by Jay Franck
c How big do you want the box to be?
c Enter 'L' (the box will go from -L to L): 
c  20
c Enter the nuclear radius (R) in units of fm: 
c  2
c Enter the number of points for the solution: 
c  200
c Data written to: 'ground.dat', 'first.dat',
c  'second.dat', and 'third.dat'
c The first four energy levels for R=    10.000000    
c Ground State Energy (MeV)=   -49.533504    
c First Excited Energy (MeV)=   -48.139069    
c Second Excited Energy (MeV)=   -45.823936    
c Third Excited Energy (MeV)=   -42.602100  





