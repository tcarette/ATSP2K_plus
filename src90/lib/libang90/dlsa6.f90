!
!     --------------------------------------------------------------
!     D L S A 6
!     --------------------------------------------------------------
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE DLSA6(K, K4, K3, K5, K2, K1, J12, IRE, IAT, REC) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:15:35  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ixjtik_I 
      USE sixj_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: K 
      INTEGER  :: K4 
      INTEGER  :: K3 
      INTEGER  :: K5 
      INTEGER  :: K2 
      INTEGER  :: K1 
      INTEGER  :: J12 
      INTEGER , INTENT(IN) :: IRE 
      INTEGER , INTENT(OUT) :: IAT 
      REAL(DOUBLE)  :: REC 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      REC = ZERO 
      IAT = 1 
      IF (IRE == 0) THEN 
         IF (IXJTIK(K4,K2,J12,K1,K5,K3) == 0) IAT = 0 
      ELSE 
         CALL SIXJ (K4, K2, J12, K1, K5, K3, 0, REC) 
         REC = REC*DSQRT(DBLE((J12 + 1)*(K3 + 1))) 
         IF (MOD(2*K5 + K2 + K4 - J12,4) /= 0) REC = -REC 
      ENDIF 
      RETURN  
      END SUBROUTINE DLSA6 
