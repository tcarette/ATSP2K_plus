!                                                                  *
!     -------------------------------------------------------------
!      T W O 3 3 A
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
!                                              N'2 = N2 + 1        *
!                                                                  *
!                                                                  *
!     CASES 2313   + + - -        TRANSFORM TO  2133   + - + -     *
!           3231                                2133               *
!                                                                  *
!                                                                  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE TWO33A(KL1, KS1, KL2, KS2, K, IA, IB, IC, INE1, INE2, CA1, CB1&
         ) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE TRK_C 
      USE TRK2_C 
      USE MEDEFN_C 
      USE PERMAT_C 
      USE CASEOP_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:23:13  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE a1a2w3lsp_I 
      USE rlsp3_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: KL1 
      INTEGER  :: KS1 
      INTEGER  :: KL2 
      INTEGER  :: KS2 
      INTEGER , INTENT(IN) :: K 
      INTEGER  :: IA 
      INTEGER  :: IB 
      INTEGER  :: IC 
      INTEGER  :: INE1 
      INTEGER  :: INE2 
      REAL(DOUBLE) , INTENT(OUT) :: CA1 
      REAL(DOUBLE) , INTENT(OUT) :: CB1 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: W, RES, REL, S 
!-----------------------------------------------
      CA1 = ZERO 
      CB1 = ZERO 
!
!     CASES 2313   + + - -        TRANSFORM TO  2133   + - + -
!
      CALL A1A2W3LSP (IB, IA, IC, IC, KL2, KS2, (-HALF), HALF, HALF, (-HALF), W&
         ) 
      IF (DABS(W) > EPS) THEN 
         IF (IRS(KS1+1,KS2+1) == 0) THEN 
            IRS(KS1+1,KS2+1) = 1 
            CALL RLSP3 (3, IB, IA, IC, 1, 1, 2*KS1, 2*KS2, 2*K, RES) 
            RS(KS1+1,KS2+1) = RES 
         ELSE 
            RES = RS(KS1+1,KS2+1) 
         ENDIF 
         IF (DABS(RES) > EPS) THEN 
            IF (IRL(KL1+1,KL2+1) == 0) THEN 
               IRL(KL1+1,KL2+1) = 1 
               CALL RLSP3 (2, IB, IA, IC, 2*ID2(3), 2*ID1(3), 2*KL1, 2*KL2, 2*K&
                  , REL) 
               RL(KL1+1,KL2+1) = REL 
            ELSE 
               REL = RL(KL1+1,KL2+1) 
            ENDIF 
            IF (DABS(REL) > EPS) THEN 
               S = DBLE((2*KL1 + 1)*(2*KL2 + 1)*(2*KS1 + 1)*(2*KS2 + 1)) 
               CA1 = HALF*W*RES*REL/DSQRT(S) 
            ENDIF 
         ENDIF 
      ENDIF 
!
!     CASES 3231   + + - -        TRANSFORM TO  2133   + - + -
!
      CALL A1A2W3LSP (IB, IA, IC, IC, KL1, KS1, (-HALF), HALF, HALF, (-HALF), W&
         ) 
      IF (DABS(W) < EPS) RETURN  
      IF (IRS(KS2+1,KS1+1) == 0) THEN 
         IRS(KS2+1,KS1+1) = 1 
         CALL RLSP3 (3, IB, IA, IC, 1, 1, 2*KS2, 2*KS1, 2*K, RES) 
         RS(KS2+1,KS1+1) = RES 
      ELSE 
         RES = RS(KS2+1,KS1+1) 
      ENDIF 
      IF (DABS(RES) > EPS) THEN 
         IF (IRL(KL2+1,KL1+1) == 0) THEN 
            IRL(KL2+1,KL1+1) = 1 
            CALL RLSP3 (2, IB, IA, IC, 2*ID2(3), 2*ID1(3), 2*KL2, 2*KL1, 2*K, &
               REL) 
            RL(KL2+1,KL1+1) = REL 
         ELSE 
            REL = RL(KL2+1,KL1+1) 
         ENDIF 
         IF (DABS(REL) > EPS) THEN 
            S = DBLE((2*KL1 + 1)*(2*KL2 + 1)*(2*KS1 + 1)*(2*KS2 + 1)) 
            CB1 = HALF*W*RES*REL/DSQRT(S) 
            IF (MOD(KL1 + KL2 + KS1 + KS2,2) /= 0) CB1 = -CB1 
         ENDIF 
      ENDIF 
      RETURN  
      END SUBROUTINE TWO33A 
