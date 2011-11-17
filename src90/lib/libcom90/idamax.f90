      integer function idamax (n, dx, incx) 
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
      integer :: i, ix 
      real(double) :: dmax 
!-----------------------------------------------
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!
!
      idamax = 0 
      if (n < 1) return  
      idamax = 1 
      if (n == 1) return  
      if (incx /= 1) then 
!
!        code for increment not equal to 1
!
         ix = 1 
         dmax = dabs(dx(1)) 
         ix = ix + incx 
         do i = 2, n 
            if (dabs(dx(ix)) > dmax) then 
               idamax = i 
               dmax = dabs(dx(ix)) 
            endif 
            ix = ix + incx 
         end do 
         return  
      endif 
!
!        code for increment equal to 1
!
      dmax = dabs(dx(1)) 
      do i = 2, n 
         if (dabs(dx(i)) <= dmax) cycle  
         idamax = i 
         dmax = dabs(dx(i)) 
      end do 
      return  
      end function idamax 
