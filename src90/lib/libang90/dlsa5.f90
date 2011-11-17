!
!     -------------------------------------------------------------
!      D L S A 5
!     -------------------------------------------------------------
!                                                                  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE DLSA5(K, JA1, KA, IRE, IAT, REC) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE MEDEFN_C 
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
      INTEGER , INTENT(IN) :: K 
      INTEGER , INTENT(IN) :: JA1 
      INTEGER  :: KA 
      INTEGER , INTENT(IN) :: IRE 
      INTEGER , INTENT(OUT) :: IAT 
      REAL(DOUBLE) , INTENT(OUT) :: REC 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K1, ITI1, ITI1S, ITI, ITIS, JI 
      REAL(DOUBLE) :: A3 
!-----------------------------------------------
      REC = ZERO 
      K1 = IHSH + IHSH - 1 
      ITI1 = J1QN1(K1,K) - 1 
      ITI1S = J1QN2(K1,K) - 1 
      K1 = K1 - 1 
      IF (JA1 == IHSH) THEN 
         ITI = J1QN1(IHSH,K) - 1 
         ITIS = J1QN2(IHSH,K) - 1 
         JI = J1QN1(K1,K) - 1 
      ELSE 
         JI = J1QN1(IHSH,K) - 1 
         ITI = J1QN1(K1,K) - 1 
         ITIS = J1QN2(K1,K) - 1 
      ENDIF 
      IF (IRE == 0) THEN 
         IF (IXJTIK(KA,ITIS,ITI,JI,ITI1,ITI1S) /= 0) IAT = 1 
      ELSE 
         CALL SIXJ (KA, ITIS, ITI, JI, ITI1, ITI1S, 0, A3) 
         REC = A3*DSQRT(DBLE((ITI + 1)*(ITI1S + 1))) 
         IF (MOD(KA + JI + ITIS + ITI1,4) /= 0) REC = -REC 
         IAT = 1 
         IF (JA1 == IHSH) RETURN  
         IF (MOD(ITI + ITIS - ITI1S - ITI1 + 2*JI,4) /= 0) REC = -REC 
      ENDIF 
      RETURN  
      END SUBROUTINE DLSA5 
