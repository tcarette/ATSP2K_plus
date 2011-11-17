 
      subroutine trbak3(nm, n, nv, a, m, z) 
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
      integer , intent(in) :: nm 
      integer , intent(in) :: n 
      integer , intent(in) :: nv 
      integer , intent(in) :: m 
      real(double) , intent(in) :: a(nv) 
      real(double) , intent(inout) :: z(nm,m) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k, l, ik, iz 
      real(double) :: h, s 
!-----------------------------------------------
!
!
!     this subroutine is a translation of the algol procedure trbak3,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine forms the eigenvectors of a real symmetric
!     matrix by back transforming those of the corresponding
!     symmetric tridiagonal matrix determined by  tred3.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        nv must be set to the dimension of the array parameter a
!          as declared in the calling program dimension statement.
!
!        a contains information about the orthogonal transformations
!          used in the reduction by  tred3  in its first
!          n*(n+1)/2 positions.
!
!        m is the number of eigenvectors to be back transformed.
!
!        z contains the eigenvectors to be back transformed
!          in its first m columns.
!
!     on output
!
!        z contains the transformed eigenvectors
!          in its first m columns.
!
!     note that trbak3 preserves vector euclidean norms.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (m /= 0) then 
         if (n /= 1) then 
!
            do i = 2, n 
               l = i - 1 
               iz = (i*l)/2 
               ik = iz + i 
               h = a(ik) 
               if (h == 0.0D0) cycle  
!
               do j = 1, m 
                  ik = iz 
!
                  s = sum(a(ik+1:l+ik)*z(:l,j)) 
!     .......... double division avoids possible underflow ..........
                  s = (s/h)/h 
                  ik = iz 
!
                  z(:l,j) = z(:l,j) - s*a(ik+1:l+ik) 
!
               end do 
!
            end do 
!
         endif 
      endif 
      return  
      end subroutine trbak3 
