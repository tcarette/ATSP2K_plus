      real(kind(0.0d0)) function dasum (n, dx, incx) 
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
      real(double) , intent(in) :: dx(1) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ix, m, mp1 
      real(double) :: dtemp 
!-----------------------------------------------
!
!     takes the sum of the absolute values.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified to correct problem with negative increment, 8/21/90.
!
!
      dasum = 0.0D0 
      dtemp = 0.0D0 
      if (n <= 0) return  
      if (incx /= 1) then 
!
!        code for increment not equal to 1
!
         ix = 1 
         if (incx < 0) ix = ((-n) + 1)*incx + 1 
         do i = 1, n 
            dtemp = dtemp + dabs(dx(ix)) 
            ix = ix + incx 
         end do 
         dasum = dtemp 
         return  
      endif 
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
      m = mod(n,6) 
      if (m == 0) go to 40 
      do i = 1, m 
         dtemp = dtemp + dabs(dx(i)) 
      end do 
      if (n < 6) go to 60 
   40 continue 
      mp1 = m + 1 
      do i = mp1, n, 6 
         dtemp = dtemp + dabs(dx(i)) + dabs(dx(i+1)) + dabs(dx(i+2)) + dabs(dx(&
            i+3)) + dabs(dx(i+4)) + dabs(dx(i+5)) 
      end do 
   60 continue 
      dasum = dtemp 
      return  
      end function dasum 
