!
!     ------------------------------------------------------------------
!     S E T V E C
!     ------------------------------------------------------------------
!
      SUBROUTINE SETVEC(VECTOR, VALUE, NDIM) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!
! VECTOR (*) = VALUE
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:43:44  11/20/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NDIM 
      REAL(DOUBLE) , INTENT(IN) :: VALUE 
      REAL(DOUBLE) , INTENT(OUT) :: VECTOR(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
!-----------------------------------------------
!
      VECTOR(:NDIM) = VALUE 
!
      RETURN  
      END SUBROUTINE SETVEC 
