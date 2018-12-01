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
	double precision y(0:100),x(0:100),h,k1,k2,k3,k4,dumy,dumx
	integer N,i,j


c.......OUTPUT........................................................

c.......INTERMEDIATE..................................................


cc......USER.INPUT.SECTION............................................

	print *,"This is computational project 5"
	print *,"Written by Jay Franck"
	
	do i=0,100
		y(i)=0.0
		x(i)=0.0
	enddo


cc......MAIN.PROGRAM.SECTION...........................................
	h=1.5
	N=2
	y(0)=5.0

	do i=0,N-1
		dumy=y(i)
		dumx=x(i)

		call fun(dumy,dumx,k1)
		call fun((dumy+(0.5*k1*h)),(dumx+(0.5*h)),k2)
		call fun((dumy+(0.5*k2*h)),(dumx+(0.5*h)),k3)
		call fun((dumy+(k3*h)),(dumx+h),k4)

		y(i+1)= y(i)+(1.0/6.0)*(k1+(2*k2)+(2*k3)+k4)*h

		x(i+1)=x(i)+h
	enddo
	


c.......OUTPUT.........................................................



      	return
35     	end



cc......SAMPLE OUTPUT.............



cc......SUBROUTINE................

	subroutine fun(y,x,k)

	double precision y,x,k
	
	k= 3*exp(-x) -0.4*y
c	print *,k
	return
	end


