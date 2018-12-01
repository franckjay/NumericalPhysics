	      subroutine simpson(Np,N,array,h,ans)
C
C  subroutine to compute integral using simpson's rule
C  CWJ - SDSU - August 2006
C
C INPUT:
C   Np: declared dimension of array
C   N : used dimension of array (note, must be even!) 
C   array(0:Np): array of points of function to be integrated
C   h :  dx 
C
C OUTPUT:
C   ans : integral
C
 
      implicit none

C..........INPUT...............
      integer Np, N        ! dimensions
      double precision array(0:Np)     ! this is legal, since being fed from outside
      double precision h               ! = dx
C.........OUTPUT................
      double precision ans
C.........INTERMEDIATE..........
      integer i

C............. ERROR TRAP...........
C..... check that N is even.........

      if( (n/2)*2 .ne. n)then
        print*,' Dimension of array is not even ',N
        stop
      endif
C................. END ERROR TRAP

      ans = 0.0            ! initialize output

      do i = 1,N-1,2 
        ans = ans+array(i-1)+ 4*array(i)+array(i+1)
      enddo
      ans = ans*h/3.

      return
      end
	  
	subroutine trapezoid(Np,N,array,h,ans)
C
C  subroutine to compute integral using trapezoidal rule
C  CWJ - SDSU - August 2006
C
C INPUT:
C   Np: declared dimension of array
C   N : used dimension of array
C   array(0:Np): array of points of function to be integrated
C   h :  dx 
C
C OUTPUT:
C   ans : integral
C
 
      implicit none

C..........INPUT...............
      integer Np, N        ! dimensions
      double precision array(0:Np)     ! this is legal, since being fed from outside
      double precision h               ! = dx
C.........OUTPUT................
      double precision ans
C.........INTERMEDIATE..........
      integer i


      ans = 0.0            ! initialize output
      ans = 0.5*(array(0) + array(N))
      do i = 1, N-1 
        ans = ans+array(i)
      enddo
      ans = ans*h
	  
      return
      end
