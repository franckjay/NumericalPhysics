      	subroutine boole(f,h,lp,output,maxsize)

cc	

      	implicit none
c.......TRANSFERRED.VARIABLES........................................       	
	integer maxsize	
	double precision f(0:maxsize),h,output
c.......INTERMEDIATE.................................................
	integer i,lp	
		
	output=0

c	Runs through Boole's Rule of a given array 'f' with stepsize h
	do i=0,lp,4

		output=output+((2*h/45)*((7*f(i))+(32*f(i+1))+(12*f(i+2))
     +	 	 +(32*f(i+3))+(7*f(i+4))))


	
	enddo      
  
    	return
     	end
