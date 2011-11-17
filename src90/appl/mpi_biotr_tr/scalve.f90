!
!
!     ------------------------------------------------------------------
!     S C A L V E
!     ------------------------------------------------------------------
!
      SUBROUTINE SCALVE(VECTOR, FACTOR, NDIM) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!
! CALCULATE SCALAR(FACTOR) TIMES VECTOR
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:40:50  11/20/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NDIM 
      REAL(DOUBLE) , INTENT(IN) :: FACTOR 
      REAL(DOUBLE) , INTENT(INOUT) :: VECTOR(1) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
!-----------------------------------------------
!
      VECTOR(:NDIM) = VECTOR(:NDIM)*FACTOR 
!
      RETURN  
      END SUBROUTINE SCALVE 
