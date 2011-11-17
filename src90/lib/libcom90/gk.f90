!
!     ------------------------------------------------------------------
!              G K
!     ------------------------------------------------------------------
!                             k
!       Returns the value of G (i,j).
!
!
      REAL(KIND(0.0D0)) FUNCTION GK (I, J, K, REL) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE PARAM_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ykf_I 
      USE quads_I 
      USE grad_I 
      USE quadr_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I 
      INTEGER  :: J 
      INTEGER  :: K 
      LOGICAL  :: REL 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      CALL YKF (I, J, K, REL) 
      GK = QUADS(I,J,1) 
      IF (MASS > 0) THEN 
         IF (MASS == 1) THEN 
            IF (K == 1) GK = GK + RMASS*GRAD(I,J)**2 
         ELSE 
            GK = GK*(D1 + RMASS/D2) 
            IF (K == 1) GK = GK + Z*RMASS*QUADR(I,J,1)*QUADR(J,I,-2) 
         ENDIF 
      ENDIF 
      RETURN  
      END FUNCTION GK 
