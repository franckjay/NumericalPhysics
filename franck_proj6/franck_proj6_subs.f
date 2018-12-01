	subroutine bigmat(Big,rightBig,H,N)


	integer maxsize,i,j,N
	parameter(maxsize=1000)
	real Big(maxsize,maxsize),H(maxsize,maxsize)
	real rightBig(maxsize,maxsize)

c.......Constructs the 2N Matrix

c.......Identity.................................................
	do i=1,2*N
		do j=1,2*N
			if (j .EQ. i) then
				Big(i,j) = 1.0
				rightBig(i,j) = 1.0
			endif
		
				
		enddo
	enddo
c.......Upper.Right..............................................
	do i=1,N
		do j=1,N
			
			Big(i,N+j) = -H(i,j)
			rightBig(i,N+j) = H(i,j)
			
		
				
		enddo
	enddo
c.......Bottom.Left.............................................
	do i=1,N
		do j=1,N
			
			Big(N+i,j) = H(i,j)
			rightBig(N+i,j) = -H(i,j)
							
		enddo
	enddo


	return
	end

c***************************************************************
c***************************************************************


	subroutine crank(Big,invB,N,psi)


	integer maxsize,i,j,N
	parameter(maxsize=1000)
	real Big(maxsize,maxsize),invB(maxsize,maxsize),psi(maxsize)
	real holder,temp(maxsize)


	do i=1,2*N
		temp(i)=0.0
	enddo

c	Multiplies the RHS 2Nx2N and the Wavefunction vector
	do i=1,2*N
		holder=0.0
		do j=1,2*N
			holder= holder + (Big(i,j)*psi(j))
		enddo
		temp(i)=holder
	enddo


c	Multiplies the inverse LHS 2Nx2N and the new RHS Vector
	do i=1,2*N
		holder=0.0
		do j=1,2*N
			holder= holder + (invB(i,j)*temp(j))
		enddo
c		Fills the new Wavefunction vector at this timestep
		psi(i)=holder
	enddo




	return
	end
