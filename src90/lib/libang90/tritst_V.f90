!
!     ------------------------------------------------------------------
!
!       T R I T S T
!     ------------------------------------------------------------------
!
      DOUBLE PRECISION FUNCTION TRITST (L, M, N) 
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  20:44:41  11/14/01
!...Switches:
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  20:44:41  11/14/01
!...Switches:
      IMPLICIT NONE
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:35:42  11/20/01  
!...Switches:                     
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: L 
      INTEGER, INTENT(IN) :: M 
      INTEGER, INTENT(IN) :: N 
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
      IF (N - LMN .GE. 0) THEN 
         IF (LM - N .GE. 0) THEN 
            TRITST = 0.D0 
            RETURN  
         ENDIF 
      ENDIF 
      TRITST = 1.D0 
      RETURN  
      END FUNCTION TRITST 
