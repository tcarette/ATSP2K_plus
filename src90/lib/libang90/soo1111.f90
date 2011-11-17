!
!     -------------------------------------------------------------
!      S O O 1 1 1 1
!     -------------------------------------------------------------
!                                                                  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
!                                                                  *
!                                                                  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
!
      SUBROUTINE SOO1111(IG, KL, IA) 
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
      USE soo1_I 
      USE two11_I 
      USE savenon_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IG 
      INTEGER  :: KL 
      INTEGER  :: IA 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KLP, KLM, LA 
      REAL(DOUBLE) :: A1, C1, C2, C 
!-----------------------------------------------
      KLP = KL + 1 
      KLM = KL - 1 
      LA = IJFUL(IA) 
      IF (KLM >= 0) THEN 
         CALL SOO1 (LJ(IA), KLM, 1, KL, 0, A1) 
         IF (DABS(A1) > EPS) THEN 
            CALL TWO11 (KLM, 1, KL, 0, 1, IA, C1) 
            C1 = C1*A1 
         ELSE 
            C1 = ZERO 
         ENDIF 
         CALL SOO1 (LJ(IA), KLM, 0, KL, 1, A1) 
         IF (DABS(A1) > EPS) THEN 
            CALL TWO11 (KLM, 0, KL, 1, 1, IA, C2) 
            C2 = C2*A1 
         ELSE 
            C2 = ZERO 
         ENDIF 
         C = C1 + C2 
         IF (DABS(C) > EPS) CALL SAVENON (6, C, KLM - 1, LA, LA, LA, LA, JA, JB&
            , 0) 
      ENDIF 
      IF (KLP <= IG) THEN 
         CALL SOO1 (LJ(IA), KLP, 1, KL, 0, A1) 
         IF (DABS(A1) > EPS) THEN 
            CALL TWO11 (KLP, 1, KL, 0, 1, IA, C1) 
            C1 = C1*A1 
         ELSE 
            C1 = ZERO 
         ENDIF 
         CALL SOO1 (LJ(IA), KLP, 0, KL, 1, A1) 
         IF (DABS(A1) > EPS) THEN 
            CALL TWO11 (KLP, 0, KL, 1, 1, IA, C2) 
            C2 = C2*A1 
         ELSE 
            C2 = ZERO 
         ENDIF 
         C = C1 + C2 
         IF (DABS(C) > EPS) CALL SAVENON (6, C, KL, LA, LA, LA, LA, JA, JB, 0) 
      ENDIF 
      RETURN  
      END SUBROUTINE SOO1111 
