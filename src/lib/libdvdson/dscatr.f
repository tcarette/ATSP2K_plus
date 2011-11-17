
      SUBROUTINE dscatr(n,x,incx,index,inci,y,incy)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION x(*), y(*), index(*)
********
*       Assume incx=inci=1=incy
********

      do 10 i=1,n
         y( index(i) )= x( i )
 10   continue
      return
      end
