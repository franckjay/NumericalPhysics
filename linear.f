cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Author: Jay Franck  Email: franckjay@gmail.com
cc Course: Research Instructor: Shafter
cc
cc Program Name: Linear Surface Photometry
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc PROGRAM NOTES:
cc __________________________________________________________________
cc
cc
cc SUBROUTINES:
cc __________________________________________________________________
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	program linear
	implicit none





c.......INPUTS........................................................
	double precision arcsec(0:65),mag(0:65)
c.......OUTPUT........................................................
	
c.......INTERMEDIATE..................................................
	integer i,j,k
	double precision b,m,bestm,bestb,y,temp,res,bres

	open(unit=7,file='surfPhot2403.txt',status='old')
	do i=0,200
		read (7,*,end=10) arcsec(i), mag(i)
	enddo
10	continue

	bres=1E7
	b=15.124414543969728
	m=1.6426099538813101
c	Loop over 'b'
	do i=0,1000
		b= b+0.000000001

c		Loop over 'm' slope values
		do j=0,10
			m = m+0.000000000001
			
c			Loop over arcsecond values
			
			do k=0,65
				temp = b-m*(log(arcsec(k)))
				res = res + abs(mag(k)-temp)
			
			enddo

			if (res .LT. bres) then
				bestm = m
				bestb = b
				bres = res
			endif
		enddo
	
	enddo
	print *,'Last m= ',m,'b= ',b,'res= ',res
	print *,'m= ',bestm,'b= ',bestb,'res= ',bres
	

	return
	end

c	Fortran/Phys580> gfortran linear.f
c	Fortran/Phys580> ./a.out
c 	Last m= 1.6436200549638862      b=    15.124424543969845      res= 12074630.205757877     
c	m=   1.6426099548803101      b=    15.124414543969728      res=11.953112124338531  


c	m=    1.6426099538813101      b=    15.124414444969727      res=    11.953111726851251 
