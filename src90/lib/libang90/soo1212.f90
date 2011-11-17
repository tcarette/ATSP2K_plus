!
!     -------------------------------------------------------------
!      S O O 1 2 1 2
!     -------------------------------------------------------------
!                                                                  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
!
      SUBROUTINE SOO1212(IG, KL, IA, IB) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE MEDEFN_C 
      USE DIAGNL_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:15:20  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE sooa_I 
      USE soob_I 
      USE two12_I 
      USE savenon_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IG 
      INTEGER  :: KL 
      INTEGER  :: IA 
      INTEGER  :: IB 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: L1, L2, LB, LA, KLP, KLM 
      REAL(DOUBLE) :: AP1, A1, C1, C2, C, C3, C4 
!-----------------------------------------------
      L1 = LJ(IA) 
      L2 = LJ(IB) 
      CALL SOOA (L1, L2, L1, L2, KL, AP1) 
      IF (DABS(AP1) < EPS) RETURN  
      LB = IJFUL(IB) 
      LA = IJFUL(IA) 
      KLP = KL + 1 
      KLM = KL - 1 
      IF (KLM >= 0) THEN 
         CALL SOOB (1, L1, L1, KLM, KL, A1) 
         IF (DABS(A1) > EPS) THEN 
            CALL TWO12 (KLM, 1, KL, 0, 1, IA, IB, C1) 
            CALL TWO12 (KLM, 0, KL, 1, 1, IA, IB, C2) 
            C = C1 + TWO*C2 
            C = AP1*A1*C 
            IF (DABS(C) > EPS) CALL SAVENON (6, C, KLM - 1, LB, LA, LB, LA, JA&
               , JB, 0) 
         ENDIF 
!
         CALL SOOB (1, L2, L2, KLM, KL, A1) 
         IF (DABS(A1) > EPS) THEN 
            CALL TWO12 (KL, 0, KLM, 1, 1, IA, IB, C1) 
            CALL TWO12 (KL, 1, KLM, 0, 1, IA, IB, C2) 
            C = C1 + TWO*C2 
            C = AP1*A1*C 
            IF (DABS(C) > EPS) CALL SAVENON (6, C, KLM - 1, LA, LB, LA, LB, JA&
               , JB, 0) 
         ENDIF 
      ENDIF 
      IF (KLP <= IG) THEN 
         CALL SOOB (1, L1, L1, KLP, KL, A1) 
         IF (DABS(A1) > EPS) THEN 
            CALL TWO12 (KLP, 1, KL, 0, 1, IA, IB, C3) 
            CALL TWO12 (KLP, 0, KL, 1, 1, IA, IB, C4) 
            C = C3 + TWO*C4 
            C = AP1*A1*C 
            IF (DABS(C) > EPS) CALL SAVENON (6, C, KL, LA, LB, LA, LB, JA, JB, &
               0) 
         ENDIF 
!
         CALL SOOB (1, L2, L2, KLP, KL, A1) 
         IF (DABS(A1) > EPS) THEN 
            CALL TWO12 (KL, 0, KLP, 1, 1, IA, IB, C3) 
            CALL TWO12 (KL, 1, KLP, 0, 1, IA, IB, C4) 
            C = C3 + TWO*C4 
            C = AP1*A1*C 
            IF (DABS(C) > EPS) CALL SAVENON (6, C, KL, LB, LA, LB, LA, JA, JB, &
               0) 
         ENDIF 
      ENDIF 
      RETURN  
      END SUBROUTINE SOO1212 
