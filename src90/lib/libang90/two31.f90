!                                                                  *
!     -------------------------------------------------------------
!      T W O 3 1
!     -------------------------------------------------------------
!                                                                  *
!                                                                  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
!                                              N'2 = N2 + 1        *
!                                                                  *
!                                                                  *
!     CASES 2111   + + - -        TRANSFORM TO  1112   + - - +     *
!           1211                                1112               *
!                                                                  *
!                                                                  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE TWO31(KL1, KS1, KL2, KS2, K, IA, IB, INE1, INE2, INE3, CA1, &
         CB1) 
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
      USE ixjtik_I 
      USE rlsp2_I 
      USE sixj_I 
      USE itrexg_I 
      USE a1aw2lsp_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: KL1 
      INTEGER , INTENT(IN) :: KS1 
      INTEGER , INTENT(IN) :: KL2 
      INTEGER , INTENT(IN) :: KS2 
      INTEGER  :: K 
      INTEGER , INTENT(IN) :: IA 
      INTEGER , INTENT(IN) :: IB 
      INTEGER  :: INE1 
      INTEGER  :: INE2 
      INTEGER  :: INE3 
      REAL(DOUBLE) , INTENT(OUT) :: CA1 
      REAL(DOUBLE) , INTENT(OUT) :: CB1 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IAA, IBB, IPS2, K2KS3, ISES1, ISES2, KIS1, KIS2, KIIS1, IAT, &
         I2, KKS2, ISES3, ISES4, IP1, IKK, IG1, I5, KKL3, ISES1L, ISES2L, KIL1&
         , KIL2, IP2, IG2, I4, KKL2, ISES3L, ISES4L 
      REAL(DOUBLE) :: REC, RECS, SN1, SN2, SN3, SN4, RECL, BKKS2, W, SLA1, S3, &
         SLB1, S4, CA, CB 
!-----------------------------------------------
      CA1 = ZERO 
      CB1 = ZERO 
      IAA = MIN0(IA,IB) 
      IBB = MAX0(IA,IB) 
!
      IF (IOCASE == 1) THEN 
! Spin-spin operator
         IPS2 = 2 
      ELSE 
         IPS2 = 1 
      ENDIF 
      DO K2KS3 = 1, 3, 2 
         ISES1 = IXJTIK(2*KS1,2*KS2,2*K,K2KS3,1,1) 
         ISES2 = IXJTIK(2*KS2,2*KS1,2*K,K2KS3,1,1) 
         IF (ISES1 + ISES2 == 0) CYCLE  
         IF (IA == IAA) THEN 
            KIS1 = K2KS3 
            KIS2 = 1 
         ELSE 
            KIS1 = 1 
            KIS2 = K2KS3 
         ENDIF 
         KIIS1 = 1 
         IF (K2KS3 == 3) KIIS1 = 2 
         IF (IRS(KIIS1,1) == 0) THEN 
            IRS(KIIS1,1) = 1 
            CALL RLSP2 (3, IAA, IBB, KIS1, KIS2, 2*K, 0, IAT, REC) 
            IF (IAT == 0) THEN 
               RECS = ZERO 
            ELSE 
               CALL RLSP2 (3, IAA, IBB, KIS1, KIS2, 2*K, 1, IAT, RECS) 
            ENDIF 
            RS(KIIS1,1) = RECS 
         ELSE 
            RECS = RS(KIIS1,1) 
         ENDIF 
         IF (DABS(RECS) < EPS) CYCLE  
         RECS = RECS*DSQRT(DBLE(K2KS3 + 1)) 
         IF (ISES1 /= 0) THEN 
            CALL SIXJ (2*KS1, 2*KS2, 2*K, K2KS3, 1, 1, 0, SN1) 
         ELSE 
            SN1 = ZERO 
         ENDIF 
         IF (ISES2 /= 0) THEN 
            CALL SIXJ (2*KS2, 2*KS1, 2*K, K2KS3, 1, 1, 0, SN2) 
         ELSE 
            SN2 = ZERO 
         ENDIF 
         DO I2 = IPS2, 2 
            KKS2 = I2 - 1 
            ISES3 = IXJTIK(2*KKS2,1,1,2*KS2,1,K2KS3) 
            ISES4 = IXJTIK(2*KKS2,1,1,2*KS1,1,K2KS3) 
            IF (ISES3 + ISES4 == 0) CYCLE  
            IF (ISES1*ISES3 /= 0) THEN 
               CALL SIXJ (2*KKS2, 1, 1, 2*KS2, 1, K2KS3, 0, SN3) 
               SN3 = SN3*DSQRT(DBLE(2*KKS2 + 1)) 
               IF (MOD(KKS2 + K2KS3,2) /= 0) SN3 = -SN3 
            ELSE 
               SN3 = ZERO 
            ENDIF 
            IF (ISES2*ISES4 /= 0) THEN 
               CALL SIXJ (2*KKS2, 1, 1, 2*KS1, 1, K2KS3, 0, SN4) 
               SN4 = SN4*DSQRT(DBLE(2*KKS2 + 1)) 
               IF (MOD(2*KKS2 + KS1 + KS2 + K2KS3,2) /= 0) SN4 = -SN4 
            ELSE 
               SN4 = ZERO 
            ENDIF 
            IP1 = ITREXG(K,LJ(IB),IK1(5)/2,ID1(5)/2,IKK) + 1 
            IF (IKK <= 0) CYCLE  
            IG1 = IP1 + IKK - 1 
            DO I5 = IP1, IG1 
               KKL3 = I5 - 1 
               ISES1L = IXJTIK(2*KL1,2*KL2,2*K,2*KKL3,2*ID2(3),2*ID1(3)) 
               ISES2L = IXJTIK(2*KL2,2*KL1,2*K,2*KKL3,2*ID2(3),2*ID1(3)) 
               IF (ISES1L + ISES2L == 0) CYCLE  
               IF (IA == IAA) THEN 
                  KIL1 = 2*KKL3 
                  KIL2 = 2*ID2(3) 
               ELSE 
                  KIL1 = 2*ID2(3) 
                  KIL2 = 2*KKL3 
               ENDIF 
               IF (IRL(I5,1) == 0) THEN 
                  CALL RLSP2 (2, IAA, IBB, KIL1, KIL2, 2*K, 0, IAT, REC) 
                  IF (IAT == 0) THEN 
                     RECL = ZERO 
                     RL(I5,1) = RECL 
                     IRL(I5,1) = 1 
                  ELSE 
                     RECL = ONE 
                  ENDIF 
               ELSE 
                  RECL = RL(I5,1) 
               ENDIF 
               IF (DABS(RECL) < EPS) CYCLE  
               IF (IRL(I5,1) == 0) THEN 
                  IRL(I5,1) = 1 
                  CALL RLSP2 (2, IAA, IBB, KIL1, KIL2, 2*K, 1, IAT, RECL) 
                  RL(I5,1) = RECL 
               ENDIF 
               IF (DABS(RECL) < EPS) CYCLE  
               IP2 = ITREXG(LJ(IA),LJ(IA),LJ(IA),KKL3,IKK) + 1 
               IF (IKK <= 0) CYCLE  
               IG2 = IP2 + IKK - 1 
               DO I4 = IP2, IG2 
                  KKL2 = I4 - 1 
                  IAT = 1 
                  IF (IOCASE == 1) THEN 
! Spin-spin operator
                     IF (MOD(1 + KKL2,2) /= 0) IAT = 0 
                  ENDIF 
                  IF (IAT == 0) CYCLE  
                  ISES3L = IXJTIK(2*KKL2,2*ID1(3),2*ID1(3),2*KL2,2*ID1(3),2*&
                     KKL3) 
                  ISES4L = IXJTIK(2*KKL2,2*ID1(3),2*ID1(3),2*KL1,2*ID1(3),2*&
                     KKL3) 
                  IF (ISES3L + ISES4L == 0) CYCLE  
                  BKKS2 = HALF*DBLE(K2KS3) 
                  CALL A1AW2LSP (IAA, IBB, KKL2, KKS2, KKL3, BKKS2, HALF, HALF&
                     , (-HALF), (-HALF), W) 
                  IF (DABS(W) < EPS) CYCLE  
                  IF (ISES1*ISES3*ISES1L*ISES3L /= 0) THEN 
                     CALL SIXJ (2*KL1, 2*KL2, 2*K, 2*KKL3, 2*ID2(3), 2*ID1(3), &
                        0, SLA1) 
                     CALL SIXJ (2*KKL2, 2*ID1(3), 2*ID1(3), 2*KL2, 2*ID1(3), 2*&
                        KKL3, 0, S3) 
                     SLA1 = SLA1*S3*DSQRT(DBLE(2*KKL2 + 1)) 
                     IF (MOD(KKL2 + ID1(3)+ID2(3)+2*KKL3,2) /= 0) SLA1 = -SLA1 
                  ELSE 
                     SLA1 = ZERO 
                  ENDIF 
                  IF (ISES2*ISES4*ISES2L*ISES4L /= 0) THEN 
                     CALL SIXJ (2*KL2, 2*KL1, 2*K, 2*KKL3, 2*ID2(3), 2*ID1(3), &
                        0, SLB1) 
                     CALL SIXJ (2*KKL2, 2*ID1(3), 2*ID1(3), 2*KL1, 2*ID1(3), 2*&
                        KKL3, 0, S4) 
                     SLB1 = SLB1*S4*DSQRT(DBLE(2*KKL2 + 1)) 
                     IF (MOD(2*KKL2 + KL1 + KL2 + 2*KKL3 + ID1(3)+ID2(3),2) /= &
                        0) SLB1 = -SLB1 
                  ELSE 
                     SLB1 = ZERO 
                  ENDIF 
                  CA = DSQRT(DBLE(2*KKL3 + 1))*SN1*SN3*RECS*SLA1*W*RECL 
                  CB = DSQRT(DBLE(2*KKL3 + 1))*SN2*SN4*RECS*SLB1*W*RECL 
                  IF (IA == IAA) THEN 
                     IF (MOD(KIL1 + KIS1 + KIL2 + KIS2 - 2*K - 2*K + 2,4) /= 0&
                        ) THEN 
                        CA = -CA 
                        CB = -CB 
                     ENDIF 
                  ENDIF 
                  CA1 = CA1 + CA 
                  CB1 = CB1 + CB 
               END DO 
            END DO 
         END DO 
      END DO 
      CA1 = CA1*HALF 
      CB1 = CB1*HALF 
      RETURN  
      END SUBROUTINE TWO31 
