      subroutine daxpy(n, da, dx, incx, dy, incy) 
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
      real(double) , intent(in) :: da 
      real(double) , intent(in) :: dx(1) 
      real(double) , intent(inout) :: dy(1) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ix, iy, m, mp1 
!-----------------------------------------------
!
!     constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
!
      if (n <= 0) return  
      if (da == 0.0D0) return  
      if (incx/=1 .or. incy/=1) then 
!
!        code for unequal increments or equal increments
!          not equal to 1
!
         ix = 1 
         iy = 1 
         if (incx < 0) ix = ((-n) + 1)*incx + 1 
         if (incy < 0) iy = ((-n) + 1)*incy + 1 
         dy(iy:(n-1)*incy+iy:incy) = dy(iy:(n-1)*incy+iy:incy) + da*dx(ix:(n-1)&
            *incx+ix:incx) 
         return  
      endif 
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
      m = mod(n,4) 
      if (m /= 0) then 
         dy(:m) = dy(:m) + da*dx(:m) 
         if (n < 4) return  
      endif 
      mp1 = m + 1 
      dy(mp1:((n-mp1+4)/4)*4-1+mp1) = dy(mp1:((n-mp1+4)/4)*4-1+mp1) + da*dx(mp1&
         :((n-mp1+4)/4)*4-1+mp1) 
      return  
      end subroutine daxpy 
