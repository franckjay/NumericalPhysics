c**********************************
c	program numint
c
c	CTK, 9/2011
c
c	Uses Simpson's rule to calculate the integral of 4/(1+x**2) on the
c	interval [0,1].  The analytic value is pi.
c
c	subroutine called:  trapezoid, simpson (found in quads.f)
c
c	function called: f(x) : the integrand
c	
c************************************

	program numint
		implicit none
		double precision f   !the integrand
		double precision xmin, xmax
		parameter(xmin = 0.) !lower bound of integral
		parameter(xmax = 1.) !upper bound of integral

		integer N, i, j, nsteps
		parameter (nsteps = 12)	!number of points to plot in error plot

		integer maxsize
		parameter (maxsize = 100000) !max number of points in integral

		double precision farray(0:maxsize)
		
		double precision dx
		double precision simp, trap
		
		double precision pi
		

c	accurate means to determine pi, to double precision
		pi = 2. * dasin(1.d0)
		
		
		print*, 'Theoretical value of integral of 4/(1+x**2) from ',
     $		'0 to 1 is pi = ', pi
		print*, 'Errors for trapezoidal and Simpsons rules written to ',
     $		'numint.txt.'
		
		open(unit=1, file='numint.txt')

		write(1,*) '# created by quadrun.f'
		write(1,*) '# Error in calculating integral of 4/(1+x**2) from ',
     $		'0 to 1 = pi'
		write(1,*) '# 1: dx, step size'
		write(1,*) '# 2: N, number of function evaluations'
		write(1,*) '# 3: error from trapezoidal rule'
		write(1,*) '# 4: error from Simpsons rule'
		
		do i=0, nsteps
			
			N = 2 * 2**(i)	!logarithmic spacing in N, dx

c		Error trap, N must be less than maxsize.
c		For more precision, increase maxsize
			if(N .gt. maxsize) then
				print*, 'stopping b/c N > maxsize'
				stop
			end if
			dx = (xmax - xmin) / float(N)
			
c		Fill up reg spaced lattice with function values
			do j = 0, N
				farray(j) = f(xmin + dx * float(j))
			end do

c		Use quadrature routines to determine numerical integration result
			call simpson(maxsize, N, farray, dx, simp)
			call trapezoid(maxsize, N, farray, dx, trap)
						
			write(1,*) dx, N, trap - pi, simp - pi
			
		end do
		
		
		close(1)
		
	end program
	
c	Value of f(x), the integrand
c	Input: x, the dummy variable of the integral
c	Returns: f, the value of the integrand at x.
	double precision function f(x)
	
		double precision x
		
		f = 4. / (1. + x**2)
	end

