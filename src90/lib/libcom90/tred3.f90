      subroutine tred3(n, nv, a, d, e, e2) 
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
      integer , intent(in) :: nv 
      real(double) , intent(inout) :: a(nv) 
      real(double) , intent(inout) :: d(n) 
      real(double) , intent(inout) :: e(n) 
      real(double) , intent(out) :: e2(n) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k, l, ii, iz, jk, jm1 
      real(double) :: f, g, h, hh, scale 
!-----------------------------------------------
!
!
!     this subroutine is a translation of the algol procedure tred3,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a real symmetric matrix, stored as
!     a one-dimensional array, to a symmetric tridiagonal matrix
!     using orthogonal similarity transformations.
!
!     on input
!
!        n is the order of the matrix.
!
!        nv must be set to the dimension of the array parameter a
!          as declared in the calling program dimension statement.
!
!        a contains the lower triangle of the real symmetric
!          input matrix, stored row-wise as a one-dimensional
!          array, in its first n*(n+1)/2 positions.
!
!     on output
!
!        a contains information about the orthogonal
!          transformations used in the reduction.
!
!        d contains the diagonal elements of the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2 may coincide with e if the squares are not needed.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
!     .......... for i=n step -1 until 1 do -- ..........
      do ii = 1, n 
         i = n + 1 - ii 
         l = i - 1 
         iz = (i*l)/2 
         h = 0.0D0 
         scale = 0.0D0 
         if (l < 1) go to 130 
!     .......... scale row (algol tol then not needed) ..........
         do k = 1, l 
            iz = iz + 1 
            d(k) = a(iz) 
            scale = scale + dabs(d(k)) 
         end do 
!
         if (scale /= 0.0D0) go to 140 
  130    continue 
         e(i) = 0.0D0 
         e2(i) = 0.0D0 
         go to 290 
!
  140    continue 
         d(:l) = d(:l)/scale 
         h = sum(d(:l)*d(:l)) 
!
         e2(i) = scale*scale*h 
         f = d(l) 
         g = -dsign(dsqrt(h),f) 
         e(i) = scale*g 
         h = h - f*g 
         d(l) = f - g 
         a(iz) = scale*d(l) 
         if (l /= 1) then 
            jk = 1 
!
            do j = 1, l 
               f = d(j) 
               g = 0.0D0 
               jm1 = j - 1 
               if (jm1 >= 1) then 
!
                  if (jm1 > 0) then 
                     g = sum(a(jk:jm1-1+jk)*d(:jm1)) 
                     e(:jm1) = e(:jm1) + a(jk:jm1-1+jk)*f 
                     jk = jm1 + jk 
                  endif 
               endif 
!
               e(j) = g + a(jk)*f 
               jk = jk + 1 
            end do 
!     .......... form p ..........
            f = 0.0D0 
!
            e(:l) = e(:l)/h 
            f = sum(e(:l)*d(:l)) 
!
            hh = f/(h + h) 
!     .......... form q ..........
            e(:l) = e(:l) - hh*d(:l) 
!
            jk = 1 
!     .......... form reduced a ..........
            do j = 1, l 
               f = d(j) 
               g = e(j) 
!
               if (j > 0) then 
                  a(jk:j-1+jk) = a(jk:j-1+jk) - f*e(:j) - g*d(:j) 
                  jk = j + jk 
               endif 
!
            end do 
!
         endif 
  290    continue 
         d(i) = a(iz+1) 
         a(iz+1) = scale*dsqrt(h) 
      end do 
!
      return  
      end subroutine tred3 
