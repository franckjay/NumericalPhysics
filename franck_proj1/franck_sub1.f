      	subroutine threepoint (xo, h,d2fx)

cc	This subroutine is called by franck_proj1.f and uses the three	
cc	point method to calculate the second derivative of the function
cc	f(x)=x*sin(x) and return the value to the user. Inputs are the
cc	value 'xo' for which the second derivative to be evaluated at
cc	and the stepsize 'h'. Output is 'd2fx'

      	implicit none
c.......TRANSFERRED.VARIABLES........................................       	
	double precision xo, h, d2fx
c.......INTERMEDIATE.................................................	
	double precision m, n	

c	For ease of reading, two of the points in the three point formula
c	are split up.
	m = xo + h
	n = xo - h

c	The Three Point Formula
	d2fx = ((m*sin(m))+(n*sin(n))-(2*(xo*sin(xo))))/(h**2)
	
      
      	   
      	return
      	end
