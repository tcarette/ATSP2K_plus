!
!     ------------------------------------------------------------------
!     W R T M A T
!     ------------------------------------------------------------------
!
      SUBROUTINE WRTMAT(A, NROW, NCOL, NMROW, NMCOL) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:49:11  11/20/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NROW 
      INTEGER , INTENT(IN) :: NCOL 
      INTEGER , INTENT(IN) :: NMROW 
      INTEGER , INTENT(IN) :: NMCOL 
      REAL(DOUBLE) , INTENT(IN) :: A(NMROW,NMCOL) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J 
!-----------------------------------------------
!
      DO I = 1, NROW 
         WRITE (6, 1010) I, (A(I,J),J=1,NCOL) 
!mrg 1010 FORMAT(1H0,I3,2X,4(1X,E14.8),/,(1H ,5X,4(1X,E14.8)))
 1010    FORMAT('0',I5,2X,4(1X,E14.8),/,(' ',7X,4(1X,E14.8))) 
      END DO 
      RETURN  
      END SUBROUTINE WRTMAT 
