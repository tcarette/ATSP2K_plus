!
!     ------------------------------------------------------------------
!     P O L I N T
!     ------------------------------------------------------------------
!
      subroutine polint(xa, ya, n, x, y, dy) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:46:35  11/20/01  
!...Switches:                     
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      real(double) , intent(in) :: x 
      real(double) , intent(out) :: y 
      real(double) , intent(out) :: dy 
      real(double) , intent(in) :: xa(n) 
      real(double) , intent(in) :: ya(n) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: nmax = 10 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, m, ns 
      real(double) :: den, dif, dift, ho, hp, w 
      real(double), dimension(nmax) :: c, d 
!-----------------------------------------------
      ns = 1 
      dif = abs(x - xa(1)) 
      do i = 1, n 
         dift = abs(x - xa(i)) 
         if (dift < dif) then 
            ns = i 
            dif = dift 
         endif 
         c(i) = ya(i) 
         d(i) = ya(i) 
      end do 
      y = ya(ns) 
      ns = ns - 1 
      do m = 1, n - 1 
         do i = 1, n - m 
            ho = xa(i) - x 
            hp = xa(i+m) - x 
            w = c(i+1) - d(i) 
            den = ho - hp 
            if (den == 0.D0) then 
               write (*, '(2A)') 'PAUSE ', 'failure in polint' 
               read * 
            endif 
            den = w/den 
            d(i) = hp*den 
            c(i) = ho*den 
         end do 
         if (2*ns < n - m) then 
            dy = c(ns+1) 
         else 
            dy = d(ns) 
            ns = ns - 1 
         endif 
         y = y + dy 
      end do 
      return  
      end subroutine polint 
