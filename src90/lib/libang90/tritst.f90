!
!     ------------------------------------------------------------------
!
!       T R I T S T
!     ------------------------------------------------------------------
!
      REAL(KIND(0.0D0)) FUNCTION TRITST (L, M, N) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  20:44:41  11/14/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: L 
      INTEGER , INTENT(IN) :: M 
      INTEGER , INTENT(IN) :: N 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LMN, LM 
!-----------------------------------------------
!
!
!      IF  TRITST=1.0   THE TRIANGLE RELATION IS NOT SATISFIED
!      IF  TRITST=0.0   THE TRIANGLE RELATION IS SATISFIED
!
      LMN = IABS(L - M) 
      LM = L + M 
      IF (N - LMN >= 0) THEN 
         IF (LM - N >= 0) THEN 
            TRITST = 0.D0 
            RETURN  
         ENDIF 
      ENDIF 
      TRITST = 1.D0 
      RETURN  
      END FUNCTION TRITST 
