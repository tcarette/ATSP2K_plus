!----------------------------------------------------------------------
!        U K
!---------------------------------------------------------------------
 
 
      REAL(KIND(0.0D0)) FUNCTION UK (I, J, II, JJ, K) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE PARAM_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  23:08:12  11/15/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dzk_I 
      USE quads_I 
      USE ykk_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I 
      INTEGER  :: J 
      INTEGER  :: II 
      INTEGER  :: JJ 
      INTEGER  :: K 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: UK1, UK2 
!-----------------------------------------------
!
      CALL DZK (J, JJ, K) 
      UK1 = QUADS(I,II,2) 
!
      CALL YKK (J, JJ, K, 1) 
      UK2 = QUADS(I,II,2) 
!
      UK = (UK1 - (K + 2.)/(K + K + 1.)*UK2)*FINE 
      RETURN  
      END FUNCTION UK 
