!
!     ------------------------------------------------------------------
!     S L O P E
!     ------------------------------------------------------------------
!
      subroutine slope 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use radial_C
      use nel_C
      use fo_C
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:45:37  11/20/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use polint_I 
      implicit none
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, i 
      real(double), dimension(4) :: xa, ya 
      real(double) :: daz 
!-----------------------------------------------
!      PARAMETER (NWD=128,NOD=220)
!      character*3 el
!      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
!      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
!     :       (IQMAX,MAX(1))
!      COMMON /NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7),el(nwd)
!      COMMON /FO/IWF(2),NCLOS(2),NOPEN(2),MAXNFO
!     Determine the starting parameter az =  p(r)/r^(l+1)|r=0
      do j = 1, maxnfo 
         do i = 4, 1, -1 
            xa(i) = r(i) 
            ya(i) = p(i,j)*dsqrt(r(i))/dexp((l(j)+1.D0)*dlog(r(i))) 
         end do 
         call polint (xa, ya, 4, 0.D0, az(j), daz) 
      end do 
      return  
      end subroutine slope 
