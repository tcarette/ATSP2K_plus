      subroutine dcopy(n, dx, incx, dy, incy) 
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
      integer , intent(in) :: incy 
      real(double) , intent(in) :: dx(1) 
      real(double) , intent(out) :: dy(1) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ix, iy, m, mp1 
!-----------------------------------------------
!
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
!
      if (n <= 0) return  
      if (incx/=1 .or. incy/=1) then 
!
!        code for unequal increments or equal increments
!          not equal to 1
!
         ix = 1 
         iy = 1 
         if (incx < 0) ix = ((-n) + 1)*incx + 1 
         if (incy < 0) iy = ((-n) + 1)*incy + 1 
         dy(iy:(n-1)*incy+iy:incy) = dx(ix:(n-1)*incx+ix:incx) 
         return  
      endif 
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
      m = mod(n,7) 
      if (m /= 0) then 
         dy(:m) = dx(:m) 
         if (n < 7) return  
      endif 
      mp1 = m + 1 
      dy(mp1:((n-mp1+7)/7)*7-1+mp1) = dx(mp1:((n-mp1+7)/7)*7-1+mp1) 
      return  
      end subroutine dcopy 
