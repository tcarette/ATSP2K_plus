      subroutine dswap(n, dx, incx, dy, incy) 
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
      real(double) , intent(inout) :: dx(1) 
      real(double) , intent(inout) :: dy(1) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ix, iy, m, mp1 
      real(double) :: dtemp 
!-----------------------------------------------
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal one.
!     jack dongarra, linpack, 3/11/78.
!
!
      if (n <= 0) return  
      if (incx/=1 .or. incy/=1) then 
!
!       code for unequal increments or equal increments not equal
!         to 1
!
         ix = 1 
         iy = 1 
         if (incx < 0) ix = ((-n) + 1)*incx + 1 
         if (incy < 0) iy = ((-n) + 1)*incy + 1 
         do i = 1, n 
            dtemp = dx(ix) 
            dx(ix) = dy(iy) 
            dy(iy) = dtemp 
            ix = ix + incx 
            iy = iy + incy 
         end do 
         return  
      endif 
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
      m = mod(n,3) 
      if (m /= 0) then 
         do i = 1, m 
            dtemp = dx(i) 
            dx(i) = dy(i) 
            dy(i) = dtemp 
         end do 
         if (n < 3) return  
      endif 
      mp1 = m + 1 
      do i = mp1, n, 3 
         dtemp = dx(i) 
         dx(i) = dy(i) 
         dy(i) = dtemp 
         dtemp = dx(i+1) 
         dx(i+1) = dy(i+1) 
         dy(i+1) = dtemp 
         dtemp = dx(i+2) 
         dx(i+2) = dy(i+2) 
         dy(i+2) = dtemp 
      end do 
      return  
      end subroutine dswap 
