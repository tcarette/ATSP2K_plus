!
!     ------------------------------------------------------------------
!     T R P M A T
!     ------------------------------------------------------------------
!
      SUBROUTINE TRPMAT(XIN, NROW, NCOL, XOUT) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!
! XOUT(I,J) = XIN(J,I)
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:43:44  11/20/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NROW 
      INTEGER , INTENT(IN) :: NCOL 
      REAL(DOUBLE) , INTENT(IN) :: XIN(NROW,NCOL) 
      REAL(DOUBLE) , INTENT(OUT) :: XOUT(NCOL,NROW) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IROW, ICOL 
!-----------------------------------------------
!
      XOUT = TRANSPOSE(XIN) 
!
      RETURN  
      END SUBROUTINE TRPMAT 
