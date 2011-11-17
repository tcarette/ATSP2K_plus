!
!     -------------------------------------------------------------
!      T W O 4 2
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 + 1        *
!                                              N'2 = N2 + 1        *
!                                              N'3 = N3 - 2,       *
!                                                                  *
!     CASES 1233   + + - -        TRANSFORM TO  1233     + + - -   *
!           2133                                1233               *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE TWO42(KL1, KS1, KL2, KS2, K, IA, IB, IC, INE1, INE2, CA1, CB1) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE MEDEFN_C 
      USE TRK_C 
      USE TRK2_C 
      USE PERMAT_C 
      USE CASEOP_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:23:13  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE itrexg_I 
      USE niness_I 
      USE rlsp3_I 
      USE a1a2w3lsp_I 
      USE ninels_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: KL1 
      INTEGER  :: KS1 
      INTEGER , INTENT(IN) :: KL2 
      INTEGER  :: KS2 
      INTEGER  :: K 
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
      INTEGER :: IP1, IKK, IG1, IPS1, IPS2, I1, KKS1, I2, KKS2, IP, IG, I3, &
         KKL2, I4, KKL1, IAT, INA, INB 
      REAL(DOUBLE) :: SN1, RECS, W, S, RECL, SLA1, SLB1 
!-----------------------------------------------
      CA1 = ZERO 
      CB1 = ZERO 
      IF (IHSH <= 2) RETURN  
      IP1 = ITREXG(LJ(IC),LJ(IC),IK3(5)/2,ID3(5)/2,IKK) + 1 
      IF (IKK <= 0) RETURN  
      IG1 = IP1 + IKK - 1 
!
      IF (IOCASE == 1) THEN 
! Spin-spin operator
         IPS1 = 2 
         IPS2 = 2 
      ELSE 
         IPS1 = 1 
         IPS2 = 1 
      ENDIF 
      DO I1 = IPS1, 2 
         KKS1 = I1 - 1 
         DO I2 = IPS2, 2 
            KKS2 = I2 - 1 
            CALL NINESS (KS1, KS2, KKS1, KKS2, K, SN1) 
            IF (DABS(SN1) < EPS) CYCLE  
            SN1 = SN1*DSQRT(DBLE((2*KKS1 + 1)*(2*KKS2 + 1))) 
            IF (IRS(I1,I2) == 0) THEN 
               IRS(I1,I2) = 1 
               CALL RLSP3 (3, IA, IB, IC, 1, 1, 2*KKS1, 2*KKS2, 2*K, RECS) 
               RS(I1,I2) = RECS 
            ELSE 
               RECS = RS(I1,I2) 
            ENDIF 
            IF (DABS(RECS) < EPS) CYCLE  
            IP = IABS(LJ(IB)-LJ(IA)) + 1 
            IG = LJ(IB) + LJ(IA) + 1 
            DO I3 = IP1, IG1 
               KKL2 = I3 - 1 
               CALL A1A2W3LSP (IC, IA, IB, IC, KKL2, KKS2, HALF, HALF, (-HALF)&
                  , (-HALF), W) 
               IF (DABS(W) < EPS) CYCLE  
               DO I4 = IP, IG 
                  KKL1 = I4 - 1 
                  IAT = 1 
                  IF (IOCASE == 1) THEN 
! Spin-spin operator
                     IF (MOD(1 + KKL2,2) /= 0) IAT = 0 
                  ENDIF 
                  IF (IAT == 0) CYCLE  
                  CALL NINELS (2*ID1(3), 2*ID3(3), 2*KL1, 2*ID2(3), 2*ID3(3), 2&
                     *KL2, 2*KKL1, 2*KKL2, 2*K, 1, INA, S) 
                  CALL NINELS (2*ID2(3), 2*ID3(3), 2*KL1, 2*ID1(3), 2*ID3(3), 2&
                     *KL2, 2*KKL1, 2*KKL2, 2*K, 1, INB, S) 
                  IF (INA + INB == 0) CYCLE  
                  IF (IRL(I4,I3) == 0) THEN 
                     IRL(I4,I3) = 1 
                     CALL RLSP3 (2, IA, IB, IC, 2*ID1(3), 2*ID2(3), 2*KKL1, 2*&
                        KKL2, 2*K, RECL) 
                     RL(I4,I3) = RECL 
                  ELSE 
                     RECL = RL(I4,I3) 
                  ENDIF 
                  IF (DABS(RECL) < EPS) CYCLE  
                  IF (INA /= 0) THEN 
                     CALL NINELS (2*ID1(3), 2*ID3(3), 2*KL1, 2*ID2(3), 2*ID3(3)&
                        , 2*KL2, 2*KKL1, 2*KKL2, 2*K, 0, INA, SLA1) 
                     SLA1 = SLA1*DSQRT(DBLE((2*KKL1 + 1)*(2*KKL2 + 1))) 
                  ELSE 
                     SLA1 = ZERO 
                  ENDIF 
                  IF (INB /= 0) THEN 
                     CALL NINELS (2*ID2(3), 2*ID3(3), 2*KL1, 2*ID1(3), 2*ID3(3)&
                        , 2*KL2, 2*KKL1, 2*KKL2, 2*K, 0, INB, SLB1) 
                     SLB1 = SLB1*DSQRT(DBLE((2*KKL1 + 1)*(2*KKL2 + 1))) 
                     IF (MOD(ID1(3)+ID2(3)+KKL1+KKS1,2) /= 0) SLB1 = -SLB1 
                  ELSE 
                     SLB1 = ZERO 
                  ENDIF 
                  CA1 = CA1 + SN1*RECS*SLA1*W*RECL 
                  CB1 = CB1 + SN1*RECS*SLB1*W*RECL 
               END DO 
            END DO 
         END DO 
      END DO 
      CA1 = -CA1*HALF 
      CB1 = -CB1*HALF 
      RETURN  
      END SUBROUTINE TWO42 
