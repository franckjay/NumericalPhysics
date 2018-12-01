cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Author: Jay Franck  Email: franckjay@gmail.com
cc Course: Physics 580 Instructor: Chad Kishimoto
cc
cc Program Name: Project 2 - Integration
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc PROGRAM NOTES:
cc __________________________________________________________________
cc
cc This program performs the integration to find the flux of neutrons
cc radiating of a box detected at a detector located at xo,yo, outside
cc of the box using Boole's Method. The integral is performed in three
cc steps: first the dz integral, then the dy, and then dx. An array of 
cc values is passed to the Boole subroutine, which calculates the
cc integral.
cc
cc
cc
cc
cc SUBROUTINES:
cc __________________________________________________________________
cc franck_boole.f: Calls Boole's Rule for a given array 'f', stepsize
cc  	'h', and gives the result 'output'.
cc  	-call boole(f,h,output,maxsize)	 
cc
cc VARIABLES
cc __________________________________________________________________
cc INPUT.............................................................
cc W: Width of the box
cc H: Height
cc D: Depth of the box
cc xo: detector position along the x axis (parallel to D)
cc yo: detector position along the y axis (parallel to W)
cc lpx: lattice points in the x-direction
cc lpy: lattice points in the y-direction
cc lpz: lattice points in the z-direction
cc
cc OUTPUT............................................................
cc flux: Neutron flux at the detector (xo,yo)
cc approx: Approximation with xo>>D (point source)
cc exact: Exact flux found using regression at 
cc 	http://xuru.org/rt/PowR.asp
cc diff: Difference between 'flux' and exact 'answers'
cc INTERMEDIATE......................................................
cc a: String too many characters and needed to be broken up.
cc dx: Step Size of X
cc dy: Step Size of y	
cc dz: Step Size of z
cc g: Array for the first integral (dz)
cc f: Array for the second integral (dy)
cc phi: Array for the third/final integral (dx)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	program integration

	implicit none

c.......INPUTS........................................................
	double precision w, h, d, xo, yo
	integer lpx,lpy,lpz
c.......OUTPUT........................................................
	double precision flux,approx,exact,diff,exact2,point
c.......INTERMEDIATE..................................................
	character*100 a
	integer i,j,k,p
	integer maxsize
	parameter (maxsize=10000)
	double precision g(0:maxsize),phi(0:maxsize),f(0:maxsize) 
	double precision dx,dy,dz,radius,PI
cc......USER.INPUT.SECTION............................................
	print *,"This is computational project 2"
	print *,"Written by Jay Franck"
	print *,"Type the width, height, and depth of the reactor: "
	read *,w,h,d

c	Error trap for the user input of 'W','H', and'D'	
	if (w .LE. 0) then
		print *, "Reactor cannot have a negative width!"
		GO TO 35
	endif	
	if (h .LE. 0) then
		print *, "Reactor cannot have a negative height!"
		GO TO 35
	endif	
	if (d .LE. 0) then
		print *, "Reactor cannot have a negative depth!"
		GO TO 35
	endif


	print *, "Enter the x and y coordinates of the neutron detector: "
	read *,xo,yo

c	Error trap for the user input of 'xo'	
	if (xo .LE. 0) then
		print *, "Detector must be a positive distance away!"
		GO TO 35
	endif

c	Error trap for the user input of 'yo'	
	if (yo .LE. 0) then
		print *, "Detector must be a positive distance away!"
		GO TO 35
	endif

	
	print *,"Enter the number of lattice points for the integration"
	print *,"in the x, y and z direction. Must be a multiple of 4!"
	read *,lpx,lpy,lpz

c	Error trap for the user input of 'lp'	
	a="The number of lattice points is not a multiple of 4!"	
	if (MOD(lpx,4) .NE. 0) then
		print *,a
		GO TO 35
	endif
	if (MOD(lpy,4) .NE. 0) then
		print *,a
		GO TO 35
	endif
	if (MOD(lpz,4) .NE. 0) then
		print *,a
		GO TO 35
	endif



c-------Stuff for the first plot------------------------
c	d=4
c	w=5
c	h=6

c	print*, 'Dependance on lattice spacing written to ',
c     +	'latticedependance.txt.'
		
c	open(unit=1, file='latticedependance.txt')

c	write(1,*) '# created by franck_proj2.f'
c	write(1,*) '# Depth, Width, Height of box: '
c	write(1,*) '# ',d,w,h," respectively."
c	write(1,*) '# Detector at xo=2,yo=1'
c	write(1,*) '# Lattice points in the first column'
c	write(1,*) '# Difference between calculated and exact flux', 
c     +	' in the second column'
	
	
c	do p=4,184,4
c	flux=0
c	lpx=p
c	lpy=p
c	lpz=p
c-------End of First Plot Stuff--------------------------

c-------Second Plot Stuff--------------------------------
c	d=4
c	w=5
c	h=6
c	print*, 'Dependance on xo written to ',
c    +	'inversesquareBoole.txt and inversesquarePoint.txt.'
c		
c	open(unit=2, file='inversesquareBoole.txt')

c	write(2,*) '# created by franck_proj2.f'
c	write(2,*) '# Depth, Width, Height of box: '
c	write(2,*) '# ',d,w,h," respectively."
c	write(2,*) '# Detector at yo=5, but varying xo'
c	write(2,*) '# Lattice points: 96 in each direction'
c	write(2,*) '# First Column: values of xo in meters'
c	write(2,*) '# Calculated Flux in the Second Column'

c	open(unit=3, file='inversesquarePoint.txt')

c	write(2,*) '# created by franck_proj2.f'
c	write(2,*) '# Depth, Width, Height of box: '
c	write(2,*) '# ',d,w,h," respectively."
c	write(2,*) '# Detector at yo=1, but varying xo'
c	write(2,*) '# Lattice points: 96 in each direction'
c	write(2,*) '# First Column: values of xo in meters'
c	write(2,*) '# Point-Source Flux in the Second Column'
	    

c	do p=10,100010,100
c		xo=p
c		point=((w*h*d)/(4*PI))*(1/(((h/2)**2)+
c    +		 ((xo+(d/2))**2)+((yo-(w/2))**2)))

c-------End of Second Plot Stuff----------------------------




cc......MAIN.PROGRAM.SECTION...........................................

c	Computes the stepsizes in each direction
	
	dx = d/DBLE(lpx)
	dy = w/DBLE(lpy)
	dz = h/DBLE(lpz)
c	Teaching Fortran what Pi is...ridiculous	
	PI = 2. * dasin(1.d0)




c	Loop over 'x'
	do i=0, lpx		

c	Loop over 'y'
		do j=0, lpy	
			
c	Loop over 'z'
			do k=0, lpz
							
c				Computes the radius for the integral									
				radius=((1/(4*PI))*(1/((((i*dx)+xo)**2)
     +		 		 +(((j*dy)-yo)**2)+((k*dz)**2))))
				
				g(k)=radius
				
c				Ends loop if the step reaches the size
c				of the box height.
				if ((k*dz) .GT. h) then
					print *,'chilly'
					GOTO 31
				endif
				
			enddo

31			call boole(g,dz,lpz,f(j),maxsize)

c			Ends loop if the step reaches the size
c			of the box width.
			if ((j*dy) .GT. w) then
				GOTO 11
			endif			
		enddo
		
11		call boole(f,dy,lpy,phi(i),maxsize)

c		Ends loop if the step reaches the size
c		of the box depth.
		if ((i*dx) .GT. d) then
			GOTO 21
		endif		
				
	enddo
	
21	call boole(phi,dx,lpx,flux,maxsize)


c_______Second Plot Stuff_______________________________
c	write(2,*) xo,"   ",flux
c	write(3,*) xo,"   ",point
c	enddo
c_______End Second Plot Stuff____________________________

c_______First Plot Stuff__________________________________

cc......I could not find the regression code on the course website
cc......so I used an online regression code and threw in a text file
cc......to http://xuru.org/rt/PowR.asp which gave me the result
cc......y=1.11613422*x**(-0.2251684507). I used this as my exact answer.
cc......However, I found that this result blows up after about 100 lattice
cc......points, so I think there may be an error in the regression calc.
	exact= 1.11613422*(lpx)**(-0.2251684507)

cc......The result above looked a little strange, so I decided to do a
cc......estimate using a large number of lattice points(960x960x960) as
cc......a separate test. I recieved an estimate result of 
cc......exact2 = 0.37285508918626825. I have included a second plot as
cc......latticepoints2.ps that shows this dependance.
	exact2 = 0.37285508918626825 

	diff = abs(flux-exact2)
	
c	write(1,*) lpx,"   ",diff
c	enddo
c_______End First Plot Stuff_______________________________


c.......OUTPUT.............................................
	print *,"The neutron flux is: ",flux
	print *,"The exact flux using regression is: ",exact
	print *,"The estimated flux using lots of LPs is: ",exact2
	print *,"The difference between the calculated and",
     +   "estimated flux: ",diff	

	
c	approx= (w*h*d)/(4*PI*xo**2)
c	print *,approx, (flux-approx)
	
c	close(1)
c	close(2)
c	close(3)
      	return
35     	end

cc......SAMPLE OUTPUT.............

c Phys580/franck_proj2> gfortran franck_proj2.f franck_boole.f
c Phys580/franck_proj2> ./a.out
c This is computational project 2
c Written by Jay Franck
c Type the width, height, and depth of the reactor: 
c 4 4 4
c Enter the x and y coordinates of the neutron detector: 
c 1 1
c Enter the number of lattice points for the integration
c in the x, y and z direction. Must be a multiple of 4!
c 960 960 960
c The neutron flux is:   0.41583897461834474     
c The exact flux using regression is:   0.23779515922069550     
c The difference between the two is:   4.09955262178473001E-002
c The estimated flux using lots of LPs is:   0.37484344840049744 

