!
!     ------------------------------------------------------------------
!              S N
!     ------------------------------------------------------------------
!
!                                      3              k
!       Evaluates the integral of (1/r)  P (r) P (r) Z (i, j; r)  with
!                                         i     j
!   respect to r.
!
      REAL(KIND(0.0D0)) FUNCTION SN (I, J, II, JJ, K) 
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
      USE zk_I 
      USE quads_I 
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
!-----------------------------------------------
      CALL ZK (J, JJ, K) 
      SN = QUADS(I,II,3)*FINE 
      RETURN  
      END FUNCTION SN 
