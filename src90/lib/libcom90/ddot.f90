 
      real(kind(0.0d0)) function ddot (n, dx, incx, dy, incy) 
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
      real(double) , intent(in) :: dy(1) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ix, iy, m, mp1 
      real(double) :: dtemp 
!-----------------------------------------------
!
!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
!
      ddot = 0.0D0 
      dtemp = 0.0D0 
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
         dtemp = dot_product(dx(ix:(n-1)*incx+ix:incx),dy(iy:(n-1)*incy+iy:incy&
            )) 
         ddot = dtemp 
         return  
      endif 
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
      m = mod(n,5) 
      if (m == 0) go to 40 
      dtemp = dot_product(dx(:m),dy(:m)) 
      if (n < 5) go to 60 
   40 continue 
      mp1 = m + 1 
      dtemp = dtemp + dot_product(dx(mp1:((n-mp1+5)/5)*5-1+mp1),dy(mp1:((n-mp1+&
         5)/5)*5-1+mp1)) 
   60 continue 
      ddot = dtemp 
      return  
      end function ddot 
