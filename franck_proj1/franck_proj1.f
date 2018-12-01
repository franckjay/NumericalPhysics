cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Author: Jay Franck  Email: franckjay@gmail.com
cc Course: Physics 580 Instructor: Chad Kishimoto
cc
cc Program Name: Numerical Derivatives
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Program Notes: This program is used to obtain the numerical
cc derivatives of a function. In this case, the function f(x)=xsin(x).
cc It utilizes the symmetric three-point and five-point methods to 
cc take the second derivatives of this function at a specific value. 
cc User input is the specific value 'xo' as well as the step size 'h'. 
cc The main program calls two subroutines, one for each of the methods.
cc franck_sub1.f is the Three Point subroutine, while _sub2.f is the
cc Five-Point Method. They each return the value 'd2fx', which is the
cc value 'xo' evaluated at the second derivative of the function 
cc f(x)=x*sin(x), which is f''(x)=2cos(x)-xsin(x)
cc
cc SUBROUTINES:
cc __________________________________________________________________
cc
cc franck_sub1.f : Three Point Method ('call threepoint(xo, h,d2fx)') 
cc franck_sub2.f : Five Point Method ('call fivepoint(xo, h,d2fx)')
cc
cc VARIABLES
cc __________________________________________________________________
cc INPUT.............................................................
cc xo: This is the initial value of x that the function will use
cc ho: The user input of the step size 'h'
cc
cc OUTPUT............................................................
cc d2fx: The value of the second derivative of the function at the 
cc 	user input value of 'xo'.
cc actual: The value of 2cos(x)-xsin(x) evaluated at 'xo'
cc delx: The absolute value of the difference between d2fx (numerical)
cc 	and 'actual' at the value 'xo'. This is written to a text file
cc
cc INTERMEDIATE......................................................
cc a,b: Character values for ease of reading output
cc h: The step size used to calculate the plots. Based on the value of
cc 	'ho' that was input by user. 
cc i: Counting integer for the 'do' loops.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



	program numDerivative

	implicit none

c.......INPUTS........................................................
	double precision xo, ho
c.......OUTPUT........................................................
	double precision d2fx, actual, delx
c.......INTERMEDIATE..................................................
	character*100 a, b
	double precision h
	integer i
	
	print *,"This is computational project 1"
	print *,"Written by Jay Franck"
	
c	User inputs initial value xo and the stepsize
	print *, 'Enter a value of x0 and stepsize h: '
	read *, xo, ho 

c	Error trap for the user input of 'xo'	
	if (xo .LE. 0) then
		print *, "xo is less than or equal to zero!"
		GO TO 35
	endif
c	Error trap for the user input of 'ho'		
	if (ho .LE. 0) then
		print *, "The step size is less than or equal to zero!"
		GO TO 35
	endif
	
c	Repeats the input back to the user for confirmation	
	print *, "Your initial value xo is: ",xo
	print *, "Your step-size 'h' is: ",ho
	
c	Calculates the 'actual' value of f''(x)
	actual = (2*cos(xo))-(xo*sin(xo))
	print *,"The actual value of the Second Derivative at xo: ",actual
	a='**************************************************************'
	print *,a 


c	Initializes three point method subroutine, returns d2fx,
c	and prints the value to the user.
	call threepoint(xo, ho, d2fx)
	delx = abs(d2fx-actual)
	print *,"The value using the Three-Point Method is: ",d2fx
	print *,"The difference from the Actual Value: ",delx
	b='**************************************************************'
	print *,b

c	This section begins a do loop that decreases the value of the step
c	size while calling on the subroutine to calculate 'd2fx' of varying
c	step sizes for the plot 'franck_proj1.ps'. The values of 'h' and
c	the absolute value between the numerical derivation and the actual
c	derivation are written to a text file.
	h = ho
	do i=1, 20
		call threepoint(xo, h, d2fx)
		delx = abs(d2fx-actual)
		open (unit=20, file= 'threepoint.txt', status='unknown')
		write (20,*) h, delx
		h= h/2
	enddo

c	Resets the Variables.	
	d2fx = 0
	delx = 0

	
c	Initializes five point method subroutine, returns d2fx,
c	calculates the difference between the actual f''(xo),
c	and prints the values to the user.
	call fivepoint(xo, ho, d2fx)
	delx = abs(d2fx-actual)
	print *,"The value using the Five-Point Method is: ",d2fx
	print *,"The difference from the Actual Value: ",delx	
	
c	This section begins a do loop that decreases the value of the step
c	size while calling on the subroutine to calculate 'd2fx' of varying
c	step sizes for the plot 'franck_proj1.ps'. The values of 'h' and
c	the absolute value between the numerical derivation and the actual
c	derivation are written to a text file.
	h=ho
	do i=1, 20		
		call fivepoint(xo, h, d2fx)
		delx = abs(d2fx-actual)
		open (unit=22, file= 'fivepoint.txt', status='unknown')
		write (22,*) h, delx
		h= h/2
	enddo

35	end






cc......SAMPLE OUTPUT.............

cPhys580/franck_proj1> gfortran franck_proj1.f franck_sub1.f franck_sub2.f
cc Phys580/franck_proj1> ./a.out
cc This is computational project 1
cc Written by Jay Franck
cc Enter a value of x0 and stepsize h: 
cc 6 0.1
cc Your initial value xo is:    6.0000000000000000     
cc Your step-size 'h' is:   0.10000000000000001     
cc The actual value of the Second Derivative at xo:    3.5968335624942869     
cc **************************************************************                                      
cc The value using the Three-Point Method is:    3.5922379828932973     
cc The difference from the actual value:   4.59557960098955220E-003
cc **************************************************************                                      
cc The value using the Five-Point Method is:    3.5968253078673715     
cc The difference from the actual value:   8.25462691533829229E-006






