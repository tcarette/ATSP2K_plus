      SUBROUTINE dgathr(n,x,incx,index,inci,y,incy)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION x(*), y(*), index(*)
      
********
*	Assume incx=inci=1
********

      iy=1
      do 10 i=1,n
	 y(iy)= x( index(i) )
	 iy=iy+incy
 10   continue
      return
      end
