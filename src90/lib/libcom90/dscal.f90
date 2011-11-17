      subroutine dscal(n, da, dx, incx) 
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
      real(double) , intent(in) :: da 
      real(double) , intent(inout) :: dx(1) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, m, mp1, nincx 
!-----------------------------------------------
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!
!
      if (n <= 0) return  
      if (incx /= 1) then 
!
!        code for increment not equal to 1
!
         nincx = n*incx 
         dx(:nincx:incx) = da*dx(:nincx:incx) 
         return  
      endif 
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
      m = mod(n,5) 
      if (m /= 0) then 
         dx(:m) = da*dx(:m) 
         if (n < 5) return  
      endif 
      mp1 = m + 1 
      dx(mp1:((n-mp1+5)/5)*5-1+mp1) = da*dx(mp1:((n-mp1+5)/5)*5-1+mp1) 
      return  
      end subroutine dscal 
