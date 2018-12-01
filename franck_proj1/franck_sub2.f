      	subroutine fivepoint (xo, h, d2fx)

cc	This subroutine is called by franck_proj1.f and uses the five	
cc	point method to calculate the second derivative of the function
cc	f(x)=x*sin(x) and return the value to the user. Inputs are the
cc	value 'xo' for which the second derivative to be evaluated at
cc	and the stepsize 'h'. Output is 'd2fx'

      	implicit none
c.......TRANSFERRED.VARIABLES........................................      	
	double precision xo, h, d2fx
c.......INTERMEDIATE.................................................		
	double precision q, r, s , t, a, b

c	For readability, four of the points in the five point formula
c	are split up.

	q = xo+(2*h)
	r = xo+h
	s = xo-h
	t = xo-(2*h) 

c	Two parts of the numerator are split because of character limits
	a= -(q*sin(q))+(16*(r*sin(r)))
	b= (-30*(xo*sin(xo)))+(16*(s*sin(s)))-(t*sin(t))

c	Five Point Formula	
	d2fx = (a+b)/(12*(h**2))     
      	     
      	return
      	end
