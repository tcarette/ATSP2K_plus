      real(kind(0.0d0)) function pythag (a, b) 
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
      real(double) , intent(in) :: a 
      real(double) , intent(in) :: b 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(double) :: p, r, s, t, u 
!-----------------------------------------------
!
!     finds dsqrt(a**2+b**2) without overflow or destructive underflow
!
      p = dmax1(dabs(a),dabs(b)) 
      if (p /= 0.0D0) then 
         r = (dmin1(dabs(a),dabs(b))/p)**2 
         t = 4.0D0 + r 
         do while(t /= 4.0D0) 
            s = r/t 
            u = 1.0D0 + 2.0D0*s 
            p = u*p 
            r = (s/u)**2*r 
            t = 4.0D0 + r 
         end do 
         go to 20 
      endif 
   20 continue 
      pythag = p 
      return  
      end function pythag 
