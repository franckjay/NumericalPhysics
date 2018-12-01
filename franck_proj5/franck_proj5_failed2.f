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
	double precision kx1(4),ky1(4),kvx1(4),kvy1(4)
	double precision dumx,dumy,dumVx,dumVy
	integer i,j
	

	parameter(G=1.0)	

	

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
c	print *,"Enter the mass of the second planet 'm2': "
c	read *,mn(2)

c	print *,"Enter the xy-positions of the second planet 'x2,y2': "
c	read *,x(2),y(2)

c	print *,"Enter its xy-velocities, vx2 and vy2: "
c	read *,vx(2),vy(2)
c.......Time.Steps.....................................................
	print *,"Enter the time step size 'dt': "
	read *,dt

	print *,"Enter an integer number of steps 'N': "
	read *,N


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

c	r12 = sqrt(((x(1)-x(2))**2)+((y(1)-y(2))**2))

	call rad(x(1),y(1),r(1))
	E = 0.5*mn(1)*((vx(1)**2)+(vy(1)**2))
     +	 +((G*M*mn(1))/r(1))
     

	print*, 'Planet 1 orbit written to xy1orbit.txt'		
	open(unit=1, file='xy1orbit.txt')

	write(1,*) '# created by franck_proj5.f'
	write(1,*) 'Planet 1'
	write(1,*) '# X values in First column, Y values in Second'
	write(1,*) '# X      Y'
	write(1,*) x(1),"  ",y(1)

	


c	print*, 'Planet 2 orbit written to xy2orbit.txt'		
c	open(unit=2, file='xy2orbit.txt')

c	write(2,*) '# created by franck_proj5.f'
c	write(2,*) 'Planet 2'
c	write(2,*) '# X values in First column, Y values in Second'
c	write(2,*) '# X      Y'
c	write(2,*) x(2),"  ",y(2)

	print*, 'Orbital Energy written to enerOrbits.txt'		
	open(unit=3, file='enerOrbits.txt')

	write(3,*) '# created by franck_proj5.f'
	write(3,*) '# Time Written in First column, Energy in Second'
	write(3,*) '# TIME       ENERGY'
	write(3,*) t,"  ",E




	do i=1,N

		dumx=x(1)
		dumy=y(1)
		
		dumVx=vx(1)
		dumVy=vy(1)
		
		call rad(dumx,dumy,r(1))


		call fun(dumx,r(1),G,M,dumVx)
		call fun(dumy,r(1),G,M,dumVy)

		kx1(1)= dumVx
		ky1(1)= dumVy

		dumx=0.5*kx1(1)+dumx+(dumVx*(t+0.5*dt))
		dumy=0.5*ky1(1)+dumy+(dumVy*(t+0.5*dt))

		call rad(dumx+0.5*kx1(1),dumy+0.5*ky1(1),r(1))
		call fun(dumx,r(1),G,M,dumVx)
		call fun(dumy,r(1),G,M,dumVy)

		kx1(2)= dumVx
		ky1(2)= dumVy

		dumx=0.5*kx1(2)+dumx+(dumVx*(t+0.5*dt))
		dumy=0.5*ky1(2)+dumy+(dumVy*(t+0.5*dt))

		call rad(dumx+0.5*kx1(2),dumy+0.5*ky1(2),r(1))
		call fun(dumx,r(1),G,M,dumVx)
		call fun(dumy,r(1),G,M,dumVy)

		kx1(3)= dumVx
		ky1(3)= dumVy
		
		dumx=kx1(3)+dumx+(dumVx*(t+dt))
		dumy=ky1(3)+dumy+(dumVy*(t+dt))

		call rad(kx1(3)+dumx,ky1(3)+dumy,r(1))
		call fun(dumx,r(1),G,M,dumVx)
		call fun(dumy,r(1),G,M,dumVy)

		kx1(4)= dumVx
		ky1(4)= dumVy
		
		vx(1)=vx(1)+(1./6.)*(kx1(1)+2*kx1(2)+2*kx1(3)+kx1(4))*dt
		vy(1)=vy(1)+(1./6.)*(ky1(1)+2*ky1(2)+2*ky1(3)+ky1(4))*dt
		print *,"Vx: ",vx(1)
		print *,"Vy: ",vy(1)
		
		x(1)=x(1)+(vx(1)*dt)
		y(1)=y(1)+(vy(1)*dt)

		write(1,*) x(1),"  ",y(1)

		
		t=t+dt
		call rad(x(1),y(1),r(1))
		E = 0.5*mn(1)*((vx(1)**2)+(vy(1)**2))
     +	 	+((G*M*mn(1))/r(1))
		write(3,*) t,"  ",E
		

	enddo
	


c.......OUTPUT.........................................................

	close(1)
	close(2)
	close(3)


      	return
35     	end
	




cc......SAMPLE OUTPUT.............




	subroutine rad(x,y,r)

	double precision x,y,r

	r=sqrt((x**2)+(y**2))

	return
	end

cc...............................



	subroutine fun(p1,rad,G,M,dv)

	double precision p1,rad,G,M,dv

	dv = -(G*M*p1)/(rad**3)

	return
	end






