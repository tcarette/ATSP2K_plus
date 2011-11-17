      subroutine ifill(n, ialpha, ix, incx) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      integer , intent(in) :: ialpha 
      integer , intent(in) :: incx 
      integer , intent(out) :: ix(1) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, jx 
!-----------------------------------------------
!
!     fills a vector, x, with the scalar value of alpha
!
!
      if (n <= 0) return  
      if (incx == 1) then 
         ix(:n) = ialpha 
      else 
         jx = 1 
         ix(jx:(n-1)*incx+jx:incx) = ialpha 
      endif 
      return  
      end subroutine ifill 
