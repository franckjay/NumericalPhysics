cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Author: Jay Franck  Email: franckjay@gmail.com
cc Course: Physics 580 Instructor: Chad Kishimoto
cc
cc Program Name: Project 3 - Least Squares
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc PROGRAM NOTES:
cc __________________________________________________________________
cc
cc Given a data set of Temperatures and years (with errors for the
cc (temp), the best fit using various polynomials is given for the data
cc set. The values for the temperatures are put in an array, while the
cc time dependence is set up in a matrix of N x k, where N are the num
cc -ber of data points, and 'k' is the order of the least squares fit.
cc 
cc
cc SUBROUTINES:
cc __________________________________________________________________
cc
cc matinv.f: Gets the inverse matrix of an NxN matrix
cc leastsq: Is called to get the coefficients for a least square fit
cc
cc VARIABLES
cc __________________________________________________________________
cc INPUT.............................................................
cc Y(22): Recorded Temperatures for each year.
cc X(22,22): Matrix with the independent variable, the year-1900.
cc error: Associated Error on each temperature reading
cc k: Order of fit, as inputted by user.
cc
cc
cc OUTPUT............................................................
cc c(21): Stores all the coefficients for the least squares fit
cc chi: Chi^2 value for each fit
cc temp: Calculated temperature at the year 2050
cc
cc INTERMEDIATE......................................................
cc Integers i,j: Used for iteration
cc
cc
cc
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	program warming

	implicit none

c.......INPUTS........................................................

	real X(22,22),Y(22),error(22)
	integer k


c.......OUTPUT........................................................
	real c(21),chi,var
	real temp
c.......INTERMEDIATE..................................................
	integer i,j

cc......USER.INPUT.SECTION............................................

	print *,"This is computational project 3"
	print *,"Written by Jay Franck"
	
	print *,"The author's birthday is 10/15/87. Therefore, data",
     +   "set 7 will be used. Filename: set7.dat"

	print *,"Enter the order 'k' of the fit, where k=3 is",
     +   " a quadratic: "
	read *,k

c	Error trap for 'k'. Program automatically breaks if not an int
	if (k .LE. 0) then
		print *, "k must be positive!"
		GO TO 35
	endif

cc......MAIN.PROGRAM.SECTION...........................................


	open(unit=7,file='set7.dat')
c	Reads in the file, filling in the necessary values of matrix X
c	as well as the error and temperature values.		
	do i=1,22
		read (7,*,end=10) X(i,2),Y(i),error(i)
		X(i,1)=1
		X(i,2)=(X(i,2)-1900)
c		Uncomment the line below for the 'basic' project
c		error(i)=1.0		
		do j=2,k-1
			X(i,(1+j))=((X(i,2))**j)
		enddo
		
	enddo	
10	continue


c 	Calls the leastsq.f subroutine, which calculates the Coefficients	
	call leastsquares(X,C,Y,error,k,var)

c	Writes the results to the specified file
	print *,"The coefficients and data written to 'k_5.dat'"	
	open(unit=3,file='k_5.dat')
	write(3,*) '# created by franck_proj3.f'
	write(3,*) '# The order of the polynomial is k= ',k
	write(3,*) '# The coefficients are: '

	do i=1,k
		write(3,*)'# C(',i,')=',C(i)
	enddo

c	Calculates the computed chi squared values using the coefficients
	do i=1,22
		temp=0.0
		do j=1,k
			temp= temp + c(j)*(X(i,2)**(j-1))
		enddo
		chi= chi +((1/(error(i)**2))*((temp-Y(i))**2))
	enddo

c	Calculates the Reduced Chi^2
	chi = chi/(22-k)
	
	write(3,*)'# Reduced Chi-Squared of this fit is: ',chi
	
c	Writes out data points for the best fit line
	write(3,*) '# The best fit line is written below'
	write(3,*) '# First Column: Year____ Second Column: Temp'
	do i=4,98
		temp=0.0
		do j=1,k			
			
			temp= temp + c(j)*(i**(j-1))
			
		enddo
		write(3,*) (i+1900),'  ',temp,' ',sqrt(var)
	enddo

c	Calculates the temperature in 2050
	temp=0.0	
	do j=1,k		
		temp= temp + c(j)*(150**(j-1))			
	enddo
	
	

c.......OUTPUT.........................................................


	print *,"The projected temperature in 2050 is: ",temp, 'F'
	print *,'with an uncertainty (+/-): ',sqrt(var),'degrees'

	write(3,*) (2050),'  ',temp,'   ',sqrt(var)





	close(7)
	close(3)
      	return
35     	end



cc......SAMPLE OUTPUT.............


c Phys580/franck_proj3> gfortran franck_proj3.f leastsq.f matinv.f
c matinv.f:91.72:

c        IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'                       
                                                                        1
cWarning: Deleted feature: PAUSE statement at (1)
cPhys580/franck_proj3> ./a.out
c This is computational project 3
c Written by Jay Franck
c The author's birthday is 10/15/87. Therefore, dataset 7 will be used. Filename: set7.dat
c Enter the order 'k' of the fit, where k=3 is a quadratic: 
c  5
c The coefficients and data written to 'k_5.dat'
c The projected temperature in 2050 is:    76.266724     F
c with an uncertainty (+/-):   0.23762067     degrees





