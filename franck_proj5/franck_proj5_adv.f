cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Author: Jay Franck  Email: franckjay@gmail.com
cc Course: Physics 580 Instructor: Chad Kishimoto
cc
cc Program Name: Project 5 - Orbits
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc PROGRAM NOTES:
cc __________________________________________________________________
cc
cc
cc 
cc
cc SUBROUTINES:
cc __________________________________________________________________
cc
cc
cc
cc VARIABLES
cc __________________________________________________________________
cc INPUT.............................................................
cc
cc
cc
cc OUTPUT............................................................
cc
cc
cc
cc INTERMEDIATE......................................................
cc
cc
cc
cc
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	program orbits

	implicit none

c.......INPUTS........................................................

	double precision M,mn(2)
	
	
	double precision x(2),y(2),vx(2),vy(2),dt
	integer N


c.......OUTPUT........................................................

c.......INTERMEDIATE..................................................

	double precision G,h,r(2),t,E,r12
	parameter(G=1.0)
	double precision kx1(4),ky1(4),kvx1(4),kvy1(4)
	double precision kx2(4),ky2(4),kvx2(4),kvy2(4)
	double precision dx1,dy1,dVx1,dVy1
	double precision dx2,dy2,dVx2,dVy2
	integer i,j
	

		

	

cc......USER.INPUT.SECTION............................................

	print *,"This is computational project 5"
	print *,"Written by Jay Franck"

	print *,"Enter the mass of the primary object 'M': "
	read *,M

c.......First.Planet.Input...........................................
	print *,"Enter the mass of the first planet 'm1': "
	read *,mn(1)

	print *,"Enter the xy-positions of the first planet 'x1,y1': "
	read *,x(1),y(1)

	print *,"Enter its xy-velocities, vx1 and vy1: "
	read *,vx(1),vy(1)

c.......Second.Planet.Input...........................................
	print *,"Enter the mass of the second planet 'm2': "
	read *,mn(2)

	print *,"Enter the xy-positions of the second planet 'x2,y2': "
	read *,x(2),y(2)

	print *,"Enter its xy-velocities, vx2 and vy2: "
	read *,vx(2),vy(2)
c.......Time.Steps.....................................................
	print *,"Enter the time step size 'dt': "
	read *,dt

	print *,"Enter an integer number of steps 'N': "
	read *,N

cc****************************PRESET ORBITS****************************
cc................Uncomment to produce the supplied plots................

cc......TWO.Planets.SEPARATED........................
c	M=100.0
c	mn(1)=0.1
c	x(1)=1.0
c	y(1)=0.0
c	vx(1)=0.0
c	vy(1)=8.0
	
c	mn(2)=0.1
c	x(2)=30.0
c	y(2)=0.0
c	vx(2)=0.0
c	vy(2)=1.0
	
c	dt=0.001
c	N=100000

cc......Planet-Moon..orbit...........................


c	M=1.0
c	mn(1)=.1
c	x(1)=1.0
c	y(1)=0.0
c	vx(1)=0.0
c	vy(1)=1.0
	
c	mn(2)=3.7E-8
c	x(2)=1.0
c	y(2)=-0.0257
c	vx(2)=1.5
c	vy(2)=0.4
	
c	dt=0.001
c	N=10000


c	Error trap for the user inputs of 'M','m1', and'm2'
c	and the step size + number of steps.	
	if (M .LT. 0.0) then
		print *, "Primary Mass must not be negative!"
		GO TO 35
	endif

	if (mn(1) .LT. 0.0) then
		print *, "First planet's mass must not be negative!"
		GO TO 35
	endif

	if (mn(2) .LT. 0.0) then
		print *, "Second planet's Mass must not be negative!"
		GO TO 35
	endif
	if (dt .LE. 0.0) then
		print *, "Step size 'dt' must be positive!"
		GO TO 35
	endif
	if (MOD(N,1) .NE. 0.0) then
		print *, "Number of steps must be an integer!"
		GO TO 35
	endif


cc......MAIN.PROGRAM.SECTION...........................................


	t=0.0
	E=0.0


	r(1)=sqrt((x(1)**2)+(y(1)**2))
	r(2)=sqrt((x(2)**2)+(y(2)**2))
	r12 = sqrt(((x(1)-x(2))**2)+((y(1)-y(2))**2))

	E = 0.5*mn(1)*((vx(1)**2)+(vy(1)**2))+
     +	 0.5*mn(2)*((vx(2)**2)+(vy(2)**2))-((G*M*mn(1))/r(1))-
     +   ((G*M*mn(2))/r(2))-((G*M*mn(1)*mn(2))/r12)
     

	print*, 'Planet 1 orbit written to xy1orbit.txt'		
	open(unit=1, file='xy1orbit.txt')

	write(1,*) '# created by franck_proj5.f'
	write(1,*) '# Planet 1'
	write(1,*) '# M= ',M,' m1= ',mn(1),' N= ',N,'dt= ',dt
	write(1,*) '# vx1= ',vx(1),' vy1= ',vy(1)
	write(1,*) '# X values in First column, Y values in Second'
	write(1,*) '# X________________________Y'

	write(1,*) x(1),"  ",y(1)

	


	print*, 'Planet 2 orbit written to xy2orbit.txt'		
	open(unit=2, file='xy2orbit.txt')

	write(2,*) '# created by franck_proj5.f'
	write(2,*) '# Planet 2'
	write(2,*) '# M= ',M,' m2= ',mn(2),' N= ',N,' dt= ',dt
	write(2,*) '# vx2= ',vx(2),' vy2= ',vy(2)
	write(2,*) '# X values in First column, Y values in Second'
	write(2,*) '# X________________________Y'
	write(2,*) x(2),"  ",y(2)

	print*, 'Orbital Energy written to enerOrbits.txt'		
	open(unit=3, file='enerOrbits.txt')

	write(3,*) '# created by franck_proj5.f'
	write(3,*) '# Time Written in First column, Energy in Second'
	write(3,*) '#_TIME____________________________________ENERGY_'
	write(3,*) t,"  ",E




	do i=1,N

c...............k1........................		

		dx1=x(1)
		dy1=y(1)
		dx2=x(2)
		dy2=y(2)


		r(1)=sqrt((dx1**2)+(dy1**2))
		r(2)=sqrt((dx2**2)+(dy2**2))
		r12=sqrt(((dx1-dx2)**2)+((dy1-dy2)**2))
				
		call fn(dx1,dx2,mn(2),M,dVx1,r(1),r12)
		call fn(dy1,dy2,mn(2),M,dVy1,r(1),r12)
		call fn(dx2,dx1,mn(1),M,dVx2,r(2),r12)
		call fn(dy2,dy1,mn(1),M,dVy2,r(2),r12)

		kx1(1)=(vx(1))*dt
		ky1(1)=(vy(1))*dt
		kx2(1)=(vx(2))*dt
		ky2(1)=(vy(2))*dt

		kvx1(1)=dVx1*dt
		kvy1(1)=dVy1*dt
		kvx2(1)=dVx2*dt
		kvy2(1)=dVy2*dt


c...............k2........................	

		dx1=x(1)+0.5*kx1(1)
		dy1=y(1)+0.5*ky1(1)
		dx2=x(2)+0.5*kx2(1)
		dy2=y(2)+0.5*ky2(1)

		r(1)=sqrt((dx1**2)+(dy1**2))
		r(2)=sqrt((dx2**2)+(dy2**2))
		r12=sqrt(((dx1-dx2)**2)+((dy1-dy2)**2))
		
		call fn(dx1,dx2,mn(2),M,dVx1,r(1),r12)
		call fn(dy1,dy2,mn(2),M,dVy1,r(1),r12)
		call fn(dx2,dx1,mn(1),M,dVx2,r(2),r12)
		call fn(dy2,dy1,mn(1),M,dVy2,r(2),r12)


		kx1(2)=(vx(1)+0.5*kvx1(1))*dt
		ky1(2)=(vy(1)+0.5*kvy1(1))*dt
		kx2(2)=(vx(2)+0.5*kvx2(1))*dt
		ky2(2)=(vy(2)+0.5*kvy2(1))*dt

		kvx1(2)=dVx1*dt
		kvy1(2)=dVy1*dt
		kvx2(2)=dVx2*dt
		kvy2(2)=dVy2*dt

c...............k3........................
	
		dx1=x(1)+0.5*kx1(2)
		dy1=y(1)+0.5*ky1(2)
		dx2=x(2)+0.5*kx2(2)
		dy2=y(2)+0.5*ky2(2)

		r(1)=sqrt((dx1**2)+(dy1**2))
		r(2)=sqrt((dx2**2)+(dy2**2))
		r12=sqrt(((dx1-dx2)**2)+((dy1-dy2)**2))
		
		call fn(dx1,dx2,mn(2),M,dVx1,r(1),r12)
		call fn(dy1,dy2,mn(2),M,dVy1,r(1),r12)
		call fn(dx2,dx1,mn(1),M,dVx2,r(2),r12)
		call fn(dy2,dy1,mn(1),M,dVy2,r(2),r12)


		kx1(3)=(vx(1)+0.5*kvx1(2))*dt
		ky1(3)=(vy(1)+0.5*kvy1(2))*dt
		kx2(3)=(vx(2)+0.5*kvx2(2))*dt
		ky2(3)=(vy(2)+0.5*kvy2(2))*dt

		kvx1(3)=dVx1*dt
		kvy1(3)=dVy1*dt
		kvx2(3)=dVx2*dt
		kvy2(3)=dVy2*dt


c...............k4........................
	
		dx1=x(1)+kx1(3)
		dy1=y(1)+ky1(3)
		dx2=x(2)+kx2(3)
		dy2=y(2)+ky2(3)
	
		r(1)=sqrt((dx1**2)+(dy1**2))
		r(2)=sqrt((dx2**2)+(dy2**2))
		r12=sqrt(((dx1-dx2)**2)+((dy1-dy2)**2))
		       
		call fn(dx1,dx2,mn(2),M,dVx1,r(1),r12)
		call fn(dy1,dy2,mn(2),M,dVy1,r(1),r12)
		call fn(dx2,dx1,mn(1),M,dVx2,r(2),r12)
		call fn(dy2,dy1,mn(1),M,dVy2,r(2),r12)


		kx1(4)=(vx(1)+kvx1(3))*dt
		ky1(4)=(vy(1)+kvy1(3))*dt
		kx2(4)=(vx(2)+kvx2(3))*dt
		ky2(4)=(vy(2)+kvy2(3))*dt

		kvx1(4)=dVx1*dt
		kvy1(4)=dVy1*dt
		kvx2(4)=dVx2*dt
		kvy2(4)=dVy2*dt

c...............New Values Calculated.................................		
		x(1)=x(1)+(1.0/6.0)*(kx1(1)+2*kx1(2)+2*kx1(3)+kx1(4))
		y(1)=y(1)+(1.0/6.0)*(ky1(1)+2*ky1(2)+2*ky1(3)+ky1(4))
		x(2)=x(2)+(1.0/6.0)*(kx2(1)+2*kx2(2)+2*kx2(3)+kx2(4))
		y(2)=y(2)+(1.0/6.0)*(ky2(1)+2*ky2(2)+2*ky2(3)+ky2(4))
		
		vx(1)=vx(1)
     +		 +(1.0/6.0)*(kvx1(1)+2*kvx1(2)+2*kvx1(3)+kvx1(4))

		vy(1)=vy(1)
     +		 +(1.0/6.0)*(kvy1(1)+2*kvy1(2)+2*kvy1(3)+kvy1(4))

		vx(2)=vx(2)
     +		 +(1.0/6.0)*(kvx2(1)+2*kvx2(2)+2*kvx2(3)+kvx2(4))

		vy(2)=vy(2)
     +		 +(1.0/6.0)*(kvy2(1)+2*kvy2(2)+2*kvy2(3)+kvy2(4))


		write(1,*) x(1),"  ",y(1)
		write(2,*) x(2),"  ",y(2)

		
		t=t+dt

		r(1)=sqrt((x(1)**2)+(y(1)**2))
		r(2)=sqrt((x(2)**2)+(y(2)**2))
		r12 = sqrt(((x(1)-x(2))**2)+((y(1)-y(2))**2))

		E = 0.5*mn(1)*((vx(1)**2)+(vy(1)**2))+
     +	 	 0.5*mn(2)*((vx(2)**2)+(vy(2)**2))-((G*M*mn(1))/r(1))-
     +   	 ((G*M*mn(2))/r(2))-((G*M*mn(1)*mn(2))/r12)



		write(3,*) t,"  ",E
		

	enddo
	
	

c.......OUTPUT.........................................................

	close(1)
	close(2)
	close(3)


      	return
35     	end
	




cc......SAMPLE OUTPUT.............

cPhys580/franck_proj5> gfortran franck_proj5_adv.f
cPhys580/franck_proj5> ./a.out
c This is computational project 5
c Written by Jay Franck
c Enter the mass of the primary object 'M': 
c 100.0
c Enter the mass of the first planet 'm1': 
c 0.1
c Enter the xy-positions of the first planet 'x1,y1': 
c 1.0 0.0
c Enter its xy-velocities, vx1 and vy1: 
c 0.0 8.0
c Enter the mass of the second planet 'm2': 
c 0.1
c Enter the xy-positions of the second planet 'x2,y2': 
c 30.0 0.0
c Enter its xy-velocities, vx2 and vy2: 
c 0.0 1.0
c Enter the time step size 'dt': 
c 0.001
c Enter an integer number of steps 'N': 
c 100000
c Planet 1 orbit written to xy1orbit.txt
c Planet 2 orbit written to xy2orbit.txt
c Orbital Energy written to enerOrbits.txt
c Phys580/franck_proj5> 





cc...............................



	subroutine fn(p1,p2,mn,M,dv,rad,rad12)

	double precision p1,p2,rad,rad12,M,dv,mn

	dv=0.0

	dv = -((M*p1)/(rad**3))-((mn*(p1-p2))/(rad12**3))

	return
	end






