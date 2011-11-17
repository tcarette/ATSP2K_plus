!
!     -------------------------------------------------------------
!      T W O 2
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF ONE PARTICLE OPERATOR IN CASE :       N'1 = N1 - 2        *
!                                              N'2 = N2 + 2        *
!                                                                  *
!                                                     1122         *
!                                                                  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE TWO2(KL1, KS1, KL2, KS2, K, IA, IB, INE1, INE2, INE3, C1, C2) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE TRK_C 
      USE MEDEFN_C 
      USE PERMAT_C 
      USE CASEOP_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:23:13  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE itrexg_I 
      USE niness_I 
      USE rlsp2_I 
      USE ninels_I 
      USE w1w2lsp_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: KL1 
      INTEGER  :: KS1 
      INTEGER , INTENT(IN) :: KL2 
      INTEGER  :: KS2 
      INTEGER  :: K 
      INTEGER , INTENT(IN) :: IA 
      INTEGER , INTENT(IN) :: IB 
      INTEGER  :: INE1 
      INTEGER  :: INE2 
      INTEGER  :: INE3 
      REAL(DOUBLE) , INTENT(OUT) :: C1 
      REAL(DOUBLE) , INTENT(OUT) :: C2 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IP1, IKK, IG1, IP2, IG2, IIA, IIB, IPS1, IPS2, I1, KKS1, I2, &
         KKS2, KIS1, KIS2, IAT, I3, KKL1, I4, KKL2, KIL1, KIL2, IN 
      REAL(DOUBLE) :: SN1, REC, RECS, SN2, RECL, W, C 
!-----------------------------------------------
      C1 = ZERO 
      C2 = ZERO 
      IF (IHSH <= 1) RETURN  
      IF (IA == IB) RETURN  
      IP1 = ITREXG(LJ(IA),LJ(IA),IK1(5)/2,ID1(5)/2,IKK) + 1 
      IF (IKK <= 0) RETURN  
      IG1 = IP1 + IKK - 1 
      IP2 = ITREXG(LJ(IB),LJ(IB),IK2(5)/2,ID2(5)/2,IKK) + 1 
      IF (IKK <= 0) RETURN  
      IG2 = IP2 + IKK - 1 
      IIA = MIN0(IA,IB) 
      IIB = MAX0(IA,IB) 
      IF (IOCASE == 1) THEN 
! Spin - spin      operator
         IPS1 = 2 
         IPS2 = 2 
      ELSE 
! Spin-other-orbit operator
         IPS1 = 1 
         IPS2 = 1 
      ENDIF 
      DO I1 = IPS1, 2 
         KKS1 = I1 - 1 
         DO I2 = IPS2, 2 
            KKS2 = I2 - 1 
            IF (IA /= IIA) THEN 
               KIS1 = KKS2 
               KIS2 = KKS1 
            ELSE 
               KIS1 = KKS1 
               KIS2 = KKS2 
            ENDIF 
            CALL NINESS (KS1, KS2, KKS1, KKS2, K, SN1) 
            IF (DABS(SN1) < EPS) CYCLE  
            IF (IRS(KIS1+1,KIS2+1) == 0) THEN 
               IRS(KIS1+1,KIS2+1) = 1 
               CALL RLSP2 (3, IIA, IIB, 2*KIS1, 2*KIS2, 2*K, 0, IAT, REC) 
               IF (IAT == 0) THEN 
                  RS(KIS1+1,KIS2+1) = ZERO 
                  RECS = ZERO 
               ELSE 
                  CALL RLSP2 (3, IIA, IIB, 2*KIS1, 2*KIS2, 2*K, 1, IAT, RECS) 
                  RS(KIS1+1,KIS2+1) = RECS 
               ENDIF 
            ELSE 
               RECS = RS(KIS1+1,KIS2+1) 
            ENDIF 
            IF (DABS(RECS) < EPS) CYCLE  
            DO I3 = IP1, IG1 
               KKL1 = I3 - 1 
               DO I4 = IP2, IG2 
                  KKL2 = I4 - 1 
                  IAT = 1 
                  IF (IOCASE == 1) THEN 
! Spin-spin operator
                     IF (MOD(KKL1 + KKL2,2) /= 0) IAT = 0 
                  ENDIF 
                  IF (IAT == 0) CYCLE  
                  IF (KL1 == KL2) THEN 
                     IF (MOD(KKL1 + KKL2 + K,2) /= 0) IAT = 0 
                  ENDIF 
                  IF (IAT == 0) CYCLE  
                  IF (IA /= IIA) THEN 
                     KIL1 = KKL2 
                     KIL2 = KKL1 
                  ELSE 
                     KIL1 = KKL1 
                     KIL2 = KKL2 
                  ENDIF 
                  CALL NINELS (2*ID1(3), 2*ID2(3), 2*KL1, 2*ID1(3), 2*ID2(3), 2&
                     *KL2, 2*KKL1, 2*KKL2, 2*K, 1, IN, SN2) 
                  IF (IN == 0) CYCLE  
                  IF (IRL(KIL1+1,KIL2+1) == 0) THEN 
                     IRL(KIL1+1,KIL2+1) = 1 
                     CALL RLSP2 (2, IIA, IIB, 2*KIL1, 2*KIL2, 2*K, 0, IAT, REC) 
                     IF (IAT == 0) THEN 
                        RECL = ZERO 
                     ELSE 
                        CALL RLSP2 (2, IIA, IIB, 2*KIL1, 2*KIL2, 2*K, 1, IAT, &
                           RECL) 
                     ENDIF 
                     RL(KIL1+1,KIL2+1) = RECL 
                  ELSE 
                     RECL = RL(KIL1+1,KIL2+1) 
                  ENDIF 
                  IF (DABS(RECL) < EPS) CYCLE  
                  CALL W1W2LSP (KKL1, KKS1, KKL2, KKS2, HALF, HALF, (-HALF), (-&
                     HALF), W) 
                  IF (DABS(W) < EPS) CYCLE  
                  W = W*DSQRT(DBLE((2*KKL1 + 1)*(2*KKL2 + 1)*(2*KKS1 + 1)*(2*&
                     KKS2 + 1))) 
                  CALL NINELS (2*ID1(3), 2*ID2(3), 2*KL1, 2*ID1(3), 2*ID2(3), 2&
                     *KL2, 2*KKL1, 2*KKL2, 2*K, 0, IN, SN2) 
                  C = -RECS*SN1*W*RECL*SN2*HALF 
                  IF (IA /= IIA) THEN 
                     IF (MOD(KKL1 + KKS1 + KKL2 + KKS2 - K - K,2) /= 0) C = -C 
                  ENDIF 
                  C1 = C1 + C 
               END DO 
            END DO 
         END DO 
      END DO 
      RETURN  
      END SUBROUTINE TWO2 
