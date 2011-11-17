      subroutine dfill(n, alpha, dx, incx) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      integer , intent(in) :: incx 
      real(double) , intent(in) :: alpha 
      real(double) , intent(out) :: dx(1) 
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
         dx(:n) = alpha 
      else 
         jx = 1 
         dx(jx:(n-1)*incx+jx:incx) = alpha 
      endif 
      return  
      end subroutine dfill 
