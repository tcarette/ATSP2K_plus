!
!     ------------------------------------------------------------------
!     C O P V E C
!     ------------------------------------------------------------------
!
      SUBROUTINE COPVEC(FROM, TO, NDIM) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:43:44  11/20/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NDIM 
      REAL(DOUBLE) , INTENT(IN) :: FROM(NDIM) 
      REAL(DOUBLE) , INTENT(OUT) :: TO(NDIM) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
!-----------------------------------------------
!
!
      TO = FROM 
!
      RETURN  
      END SUBROUTINE COPVEC 
