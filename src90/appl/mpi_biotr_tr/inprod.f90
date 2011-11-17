!
!     ------------------------------------------------------------------
!     I N P R O D
!     ------------------------------------------------------------------
!
      REAL(KIND(0.0D0)) FUNCTION INPROD (A, B, NDIM) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!      CALCULATE SCALAR PRODUCT BETWEEN TO VECTORS A,B
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:50:35  11/20/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NDIM 
      REAL(DOUBLE) , INTENT(IN) :: A(2) 
      REAL(DOUBLE) , INTENT(IN) :: B(2) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
!-----------------------------------------------
!
      INPROD = DOT_PRODUCT(A(:NDIM),B(:NDIM)) 
!
      RETURN  
      END FUNCTION INPROD 
