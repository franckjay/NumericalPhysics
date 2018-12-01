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

	double precision G,h,k(4),r(2),t,E,r12
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

	r(1) = sqrt((x(1)**2)+(y(1)**2))
	r(2) = sqrt((x(2)**2)+(y(2)**2))
	r12 = sqrt(((x(1)-x(2))**2)+((y(1)-y(2))**2))

	E = 0.5*mn(1)*((vx(1)**2)+(vy(1)**2))
     +	 + 0.5*mn(1)*((vx(2)**2)+(vy(2)**2))+((G*M*mn(1))/r(1))
     +	 -((G*M*mn(2))/r(2))-((G*mn(1)*mn(2))/r12)

	print*, 'Planet 1 orbit written to xy1orbit.txt'		
	open(unit=1, file='xy1orbit.txt')

	write(1,*) '# created by franck_proj5.f'
	write(1,*) 'Planet 1'
	write(1,*) '# X values in First column, Y values in Second'
	write(1,*) '# X      Y'
	write(1,*) x(1),"  ",y(1)

	


	print*, 'Planet 2 orbit written to xy2orbit.txt'		
	open(unit=2, file='xy2orbit.txt')

	write(2,*) '# created by franck_proj5.f'
	write(2,*) 'Planet 2'
	write(2,*) '# X values in First column, Y values in Second'
	write(2,*) '# X      Y'
	write(2,*) x(2),"  ",y(2)

	print*, 'Orbital Energy written to enerOrbits.txt'		
	open(unit=3, file='enerOrbits.txt')

	write(3,*) '# created by franck_proj5.f'
	write(3,*) '# Time Written in First column, Energy in Second'
	write(3,*) '# TIME       ENERGY'
	write(3,*) t,"  ",E




	do i=1,N


c		Vx1
		call rk4(vx(1),vx,x,x,y,k,t,dt,G,M,mn(2),1)
c		Vy1
		call rk4(vy(1),vy,y,x,y,k,t,dt,G,M,mn(2),1)
c		Vx2
		call rk4(vx(2),vy,x,x,y,k,t,dt,G,M,mn(1),2)
c		Vy2
		call rk4(vy(2),vy,y,x,y,k,t,dt,G,M,mn(1),2)

		x(1)=x(1)+vx(1)*dt
		x(2)=x(2)+vx(2)*dt
		y(1)=y(1)+vy(1)*dt
		y(2)=y(2)+vy(2)*dt


		t= t+dt
		E=0.0
		E = 0.5*mn(1)*((vx(1)**2)+(vy(1)**2))
     +	 	 + 0.5*mn(1)*((vx(2)**2)+(vy(2)**2))+((G*M*mn(1))/r(1))
     +	 	 -((G*M*mn(2))/r(2))-((G*mn(1)*mn(2))/r12)

		write(1,*) x(1),"  ",y(1)
		write(2,*) x(2),"  ",y(2)
		write(3,*) t,"  ",E



	enddo	


c.......OUTPUT.........................................................

	close(1)
	close(2)
	close(3)

      	return
35     	end



cc......SAMPLE OUTPUT.............



cc......SUBROUTINES................

	subroutine fun(p1,p2,k,r,r12,G,M,mn)

	double precision p1,p2k,r,r12,G,M,mn

	
	
	k = -((G*M/(r**3))*p1)-((G*mn/(r12**3))*(p1-p2))

	return
	end

c......................................................................


	subroutine radius(x,y,r,r12)

	double precision x(2),y(2),r(2),r12

	r(1) = sqrt((x(1)**2)+(y(1)**2))
	r(2) = sqrt((x(2)**2)+(y(2)**2))
	r12 = sqrt(((x(1)-x(2))**2)+((y(1)-y(2))**2))

	return
	end

c......................................................................

	subroutine rk4(v,Vel,s,x,y,k,t,dt,G,M,mn,num)

	double precision v,Vel(2),s(2),x(2),y(2),k(4),dt,r(2),r12
	double precision a,b,M,mn
	integer num

	call radius(x,y,r,r12)
	call fun(s(1),s(2),k(1),r(num),r12,G,M,mn)
	a= s(1)+(0.5*k(1)*Vel(1))*(t+0.5*dt)
	b= s(2)+(0.5*k(1)*Vel(2))*(t+0.5*dt)
	call fun(a,b,k(2),r(num),r12,G,M,mn)
	a= s(1)+(0.5*k(2)*Vel(1))*(t+0.5*dt)
	b= s(2)+(0.5*k(2)*Vel(2))*(t+0.5*dt)
	call fun(a,b,k(3),r(num),r12,G,M,mn)
	a= s(1)+(1.0*k(3)*Vel(1))*(t+1.0*dt)
	b= s(2)+(1.0*k(3)*Vel(2))*(t+1.0*dt)
	call fun(a,b,k(4),r(num),r12,G,M,mn)

	v = v + (1.0/6.0)*(k1+(2*k2)+(2*k3)+k4)*dt
	

	return
	end




