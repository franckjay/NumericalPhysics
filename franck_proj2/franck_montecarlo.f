cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Author: Jay Franck  Email: franckjay@gmail.com
cc Course: Physics 580 Instructor: Chad Kishimoto
cc
cc Program Name: Project 2 - Integration (Monte Carlo)
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc PROGRAM NOTES:
cc __________________________________________________________________
cc
cc	Pick 'n' randomly distributed points from [0:W][0:H][0:D]
cc	(x1,y1,z1) to (x_n,y_n,z_n). Have three variables:x,y,z
cc	
cc	Find the average value: favg= (1/n)*SUM {f(x,y,z)} over all
cc	the randomly distributed points.
cc
cc	The integral is then (W-0)(H-O)(D-0)*favg
cc
cc SUBROUTINES:
cc __________________________________________________________________
cc ranlib.f : FUNCTION ran3(iseed) generates a pseudo-random number	
cc between 0 and 1 when fed an integer seed, 'iseed'.
cc
cc VARIABLES
cc __________________________________________________________________
cc INPUT.............................................................
cc W: Width of the box
cc H: Height
cc D: Depth of the box
cc xo: detector position along the x axis (parallel to D)
cc yo: detector position along the y axis (parallel to W)
cc
cc OUTPUT............................................................
cc flux: Number of neutrons detected at the xo,yo	 
cc
cc
cc INTERMEDIATE......................................................
cc ran3: Pseudo-Random number from the subroutine
cc avg: Average Value for the Monte Carlo Integration
cc rad: Calculated Radius for each point	
cc vals: Array for values used in computing the variance
cc variance: variance of the computed values
cc sig: Short for sigma, uncertainty in the flux
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	program montecarlo

	implicit none

c.......INPUTS........................................................
	double precision w, h, d, xo, yo, cp
	integer n,iseed
	
c.......OUTPUT........................................................
	double precision flux,point
c.......INTERMEDIATE..................................................
	real ran3
	double precision PI,avg,x,y,z,rad,vals(100000000),variance,sig
	integer i,p
cc......USER.INPUT.SECTION............................................


	print *,"This is computational project 2: Monte Carlo"
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

	print *,"Enter the number of points in the Monte Carlo Int: "
	read *,n
	
	print*,'Please enter a seed value (integer): '
        read*,iseed
        if(iseed.gt.0)iseed=-iseed

c	Error trap for the user input of 'n'	
	if (n .LE. 0) then
		print *, "Must have a positive number of points!"
		GO TO 35
	endif


cc......MAIN.PROGRAM.SECTION...........................................
	
	PI=4.D0*DATAN(1.D0)
c-------Second Plot Stuff--------------------------------
c	d=4
c	w=5
c	h=6
c	print*, 'Dependance on xo written to ',
c    +	'inversesquareMonte.txt.'
		
c	open(unit=2, file='inversesquareMonte.txt')

c	write(2,*) '# created by franck_montecarlo.f'
c	write(2,*) '# Depth, Width, Height of box: '
c	write(2,*) '# ',d,w,h," respectively."
c	write(2,*) '# Detector at yo=5, but varying xo'
c	write(2,*) '# Number of Points= 1E5'
c	write(2,*) '# First Column: values of xo in meters'
c	write(2,*) '# Calculated Flux in the Second Column'
c	write(2,*) '# Calculated Uncertainty in the Third Column'     

c	do p=1,10001,1000
c		xo=p


c-------End of Second Plot Stuff
	
	avg = 0	
	do i=1,n
		
		x=ran3(iseed)*d
		y=ran3(iseed)*w
		z=ran3(iseed)*h
		rad=((1/(4*PI))*(1/(((x+xo)**2)
     +		 +((y-yo)**2)+(z**2))))	
		vals(i)=rad	
		avg=avg+rad

	enddo
	avg = (1/DBLE(n))*avg
	flux = w*h*d*avg

	do i=1,n
		variance=variance+(vals(i)**2)
	enddo
	variance = abs(((1/n)*variance)-avg**2)
	sig=w*h*d*(SQRT(variance/n))
	

	print *,"The neutron flux is: ",flux
	print *,"With an uncertainty of (+/-): ",sig
	
c	write(2,*) xo,"   ",flux,"  ",sig
c	enddo
	
c	close(2)
	return
35	end

cc......SAMPLE OUTPUT.............
c Phys580/franck_proj2> gfortran franck_montecarlo.f ranlib.f
c Phys580/franck_proj2> ./a.out
c This is computational project 2: Monte Carlo
c Written by Jay Franck
c Type the width, height, and depth of the reactor: 
c 4 4 4
c Enter the x and y coordinates of the neutron detector: 
c 2000 1
c Enter the number of points in the Monte Carlo Int: 
c 100000
c Please enter a seed value (integer): 
c 4
c The neutron flux is:   1.27069518206005830E-006
c With an uncertainty of (+/-):   4.01829098711211422E-009


