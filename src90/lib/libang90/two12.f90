!                                                                  *
!     -------------------------------------------------------------
!      T W O 1 2
!     -------------------------------------------------------------
!                                                                  *
!                                                                  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
!                                                  N'2 = N2        *
!                                                                  *
!           1212   + + - -        TRANSFORM TO  1122   + - + -     *
!           2121                                1122               *
!                                                                  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE TWO12(KL1, KS1, KL2, KS2, K, IA, IB, C1) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE TRK_C 
      USE MEDEFN_C 
      USE CASEOP_C 
      USE PERMAT_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:23:13  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rlsp2_I 
      USE w1w2lsp_I 
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
      REAL(DOUBLE) , INTENT(OUT) :: C1 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IAT 
      REAL(DOUBLE) :: RECS, RECL, W, RECLS 
!-----------------------------------------------
      C1 = ZERO 
      IF (IRS(KS1+1,KS2+1) == 0) THEN 
         IRS(KS1+1,KS2+1) = 1 
         CALL RLSP2 (3, IA, IB, 2*KS1, 2*KS2, 2*K, 0, IAT, RECS) 
         IF (IAT == 0) THEN 
            RECS = ZERO 
         ELSE 
            CALL RLSP2 (3, IA, IB, 2*KS1, 2*KS2, 2*K, 1, IAT, RECS) 
         ENDIF 
         RS(KS1+1,KS2+1) = RECS 
      ELSE 
         RECS = RS(KS1+1,KS2+1) 
      ENDIF 
      IF (DABS(RECS) < EPS) RETURN  
      IF (IRL(KL1+1,KL2+1) == 0) THEN 
         IRL(KL1+1,KL2+1) = 1 
         CALL RLSP2 (2, IA, IB, 2*KL1, 2*KL2, 2*K, 0, IAT, RECL) 
         IF (IAT == 0) THEN 
            RECL = ZERO 
         ELSE 
            CALL RLSP2 (2, IA, IB, 2*KL1, 2*KL2, 2*K, 1, IAT, RECL) 
         ENDIF 
         RL(KL1+1,KL2+1) = RECL 
      ELSE 
         RECL = RL(KL1+1,KL2+1) 
      ENDIF 
      IF (DABS(RECL) < EPS) RETURN  
      CALL W1W2LSP (KL1, KS1, KL2, KS2, HALF, (-HALF), HALF, (-HALF), W) 
      IF (DABS(W) < EPS) RETURN  
      W = W/DSQRT(DBLE((2*KL1 + 1)*(2*KL2 + 1)*(2*KS1 + 1)*(2*KS2 + 1))) 
      RECLS = RECL*RECS 
      C1 = RECLS*W*HALF 
      RETURN  
      END SUBROUTINE TWO12 
