	subroutine leastsquares(x,c,y,error,k,var)
	
	implicit none


c.......TRANSFERRED.VARIABLES........................................  
	integer k
	real X(22,22),c(k),Y(22),error(22),inv(k,k),var
c.......INTERMEDIATE.................................................
	real holder,temp(k,k),temp2(k),temp3(k,k),a,varia(k),XT(22,22)
	integer h,i,j,indx(k)
c.......SUBROUTINES...................................................
c
c	matinv.f: Gets the inverse matrix of an NxN matrix	

c	Takes the Transpose of X
	do i=1,22
		do j=1,22
		XT(i,j)=X(j,i)				
		enddo
	enddo


c	Multiplies X^T and X	
	do i = 1,k
		do j = 1,k
			holder = 0.0
			temp(i,j)=0.0
			do h = 1,22
				a=1/(error(h)**2)
				holder=holder+(a)*(XT(i,h)*X(h,j))
				
			enddo
			temp(i,j) = holder			
		enddo	
	enddo

c	Multiplies X^T and the Temperatures
	do i=1,k
		holder=0.0
		do h=1,22
			a=1/(error(h)**2)
			holder=holder+a*(XT(i,h)*Y(h))
		enddo
		temp2(i)=holder
	enddo

c	Inverts [(X^T)(X)]	
	call MatInv(temp,inv,k,k,indx)

c	Calculates variance
	do i=1,k
		do j=1,k
			varia(i)=varia(i) + inv(i,j)*X(j,i)
		enddo
		
	enddo
	do i=1,k
		
			var= var + X(1,i)*varia(i)
	enddo
	var=abs(var)		
	
c	Multiplies [(X^T)(X)]^-1 and [(X^T)(Temperature)]
c	which results in the coefficients
	do i=1,k
		holder=0.0
		do j=1,k
			holder= holder + (inv(i,j)*temp2(j))
		enddo
		c(i)=holder
	enddo



        return
        end


