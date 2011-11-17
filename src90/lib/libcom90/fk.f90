!
!     ------------------------------------------------------------------
!              F K
!     ------------------------------------------------------------------
!                             k
!       Returns the value of F (i,j)
!
!
      REAL(KIND(0.0D0)) FUNCTION FK (I, J, K, REL) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE PARAM_C, MSOO=>MASS 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ykf_I 
      USE quads_I 
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
      CALL YKF (I, I, K, REL) 
      FK = QUADS(J,J,1) 
      IF (MSOO == 2) FK = FK*(D1 + RMASS/D2) 
      IF (MSOO == 5) FK = FK*(D1 + RMASS/D2) 
      RETURN  
      END FUNCTION FK 
