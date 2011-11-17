!
!     -------------------------------------------------------------
!      T W O 5 6
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
!                                              N'2 = N2 + 1        *
!     ( IREZ = 1 )                             N'3 = N3 - 1        *
!                                              N'4 = N4 + 1        *
!                                                                  *
!     CASES 1432   + + - -        TRANSFORM TO  1234     + - - +   *
!           4123                                1234               *
!                                                                  *
!                                                                  *
!                                                                  *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 + 1        *
!                                              N'2 = N2 - 1        *
!     ( IREZ = 2 )                             N'3 = N3 + 1        *
!                                              N'4 = N4 - 1        *
!                                                                  *
!     CASES 3214   + + - -        TRANSFORM TO  1234     - + + -   *
!           2341                                1234               *
!                                                                  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE TWO56(KL1, KS1, KL2, KS2, K, IA, IB, IC, ID, IREZ, CA1, CB1) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE MEDEFN_C 
      USE TRK_C 
      USE TRK2_C 
      USE KAMPAS_C 
      USE PERMAT_C 
      USE CASEOP_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:23:13  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE itrexg_I 
      USE ixjtik_I 
      USE rlsp4b_I 
      USE sixj_I 
      USE dlsa4_I 
      USE rlsp4a_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: KL1 
      INTEGER , INTENT(IN) :: KS1 
      INTEGER , INTENT(IN) :: KL2 
      INTEGER , INTENT(IN) :: KS2 
      INTEGER  :: K 
      INTEGER  :: IA 
      INTEGER  :: IB 
      INTEGER  :: IC 
      INTEGER  :: ID 
      INTEGER , INTENT(IN) :: IREZ 
      REAL(DOUBLE) , INTENT(OUT) :: CA1 
      REAL(DOUBLE) , INTENT(OUT) :: CB1 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: JN, JI, JIS, IPS1, IKK1, IGS1, IPS2, IP1, IG1, K2KS3, ISES1, &
         ISES2, KIIS1, IAT, I1, KKS1, ISES3, ISES4, I3, KKL3, IP2, IKK2, IG2, &
         ISES1L, ISES2L, I4, KKL1, ISES3L, ISES4L 
      REAL(DOUBLE) :: CAS1, CBS1, R, RECS1, SN1, SN2, RECS3, RECS2, SN3, SN4, &
         RECL1, S3, S4, RECL3, RECL2, SLA1, SLB1 
!-----------------------------------------------
      CA1 = ZERO 
      CB1 = ZERO 
      IF (IHSH <= 3) RETURN  
      JN = IHSH + IC - 1 
      JI = J1QN1(JN,3) - 1 
      JIS = J1QN2(JN,3) - 1 
      IPS1 = ITREXG(1,2*K,JI,JIS,IKK1) 
      IF (IKK1 == 0) RETURN  
      IF (IPS1 > 3) RETURN  
      IGS1 = IPS1 + IKK1 - 1 
      IGS1 = MIN0(3,IGS1) 
      IF (IOCASE == 1) THEN 
! Spin-spin operator
         IF (IGS1 /= 3) RETURN  
         IPS2 = 2 
      ELSE 
         IPS2 = 1 
      ENDIF 
      JI = (J1QN1(JN,2)-1)/2 
      JIS = (J1QN2(JN,2)-1)/2 
      IP1 = ITREXG(LJ(ID),K,JI,JIS,IKK1) + 1 
      IF (IKK1 == 0) RETURN  
      IG1 = IP1 + IKK1 - 1 
      IF (IW1(KS1+1,KS2+8) == 0) THEN 
         IW1(KS1+1,KS2+8) = 1 
         CAS1 = ZERO 
         CBS1 = ZERO 
         DO K2KS3 = IPS1, IGS1, 2 
            ISES1 = IXJTIK(2*KS2,2*KS1,2*K,K2KS3,1,1) 
            ISES2 = IXJTIK(2*KS1,2*KS2,2*K,K2KS3,1,1) 
            IF (ISES1 + ISES2 == 0) CYCLE  
            KIIS1 = 1 
            IF (K2KS3 == 3) KIIS1 = 2 
            IF (IW1(KIIS1,2) == 0) THEN 
               IW1(KIIS1,2) = 1 
               CALL RLSP4B (3, IC, ID, K2KS3, 1, 2*K, 0, IAT, R) 
               IF (IAT == 0) THEN 
                  RECS1 = ZERO 
               ELSE 
                  CALL RLSP4B (3, IC, ID, K2KS3, 1, 2*K, 1, IAT, RECS1) 
               ENDIF 
               RW1(KIIS1,2) = RECS1 
            ELSE 
               RECS1 = RW1(KIIS1,2) 
            ENDIF 
            IF (DABS(RECS1) < EPS) CYCLE  
            IF (ISES1 /= 0) CALL SIXJ (2*KS2, 2*KS1, 2*K, K2KS3, 1, 1, 0, SN1) 
            IF (ISES2 /= 0) CALL SIXJ (2*KS1, 2*KS2, 2*K, K2KS3, 1, 1, 0, SN2) 
            DO I1 = IPS2, 2 
               KKS1 = I1 - 1 
               ISES3 = IXJTIK(2*KKS1,1,1,2*KS1,1,K2KS3) 
               ISES4 = IXJTIK(2*KKS1,1,1,2*KS2,1,K2KS3) 
               IF (ISES3 + ISES4 == 0) CYCLE  
               IF (IRS(I1,KIIS1) == 0) THEN 
                  IRS(I1,KIIS1) = 1 
                  IAT = 0 
                  CALL DLSA4 (3, IB, IC, 2*KKS1, 1, K2KS3, 0, IAT, R) 
                  IF (IAT == 0) THEN 
                     RECS3 = ZERO 
                  ELSE 
                     CALL DLSA4 (3, IB, IC, 2*KKS1, 1, K2KS3, 1, IAT, RECS3) 
                  ENDIF 
                  RS(I1,KIIS1) = RECS3 
               ELSE 
                  RECS3 = RS(I1,KIIS1) 
               ENDIF 
               IF (DABS(RECS3) <= EPS) CYCLE  
               IF (IW1(I1,3) == 0) THEN 
                  IW1(I1,3) = 1 
                  CALL RLSP4A (3, IA, IB, IC, 1, 1, 2*KKS1, 0, IAT, R) 
                  IF (IAT == 0) THEN 
                     RECS2 = ZERO 
                  ELSE 
                     CALL RLSP4A (3, IA, IB, IC, 1, 1, 2*KKS1, 1, IAT, RECS2) 
                  ENDIF 
                  RW1(I1,3) = RECS2 
               ELSE 
                  RECS2 = RW1(I1,3) 
               ENDIF 
               IF (DABS(RECS2) <= EPS) CYCLE  
               IF (ISES1*ISES3 /= 0) THEN 
                  CALL SIXJ (2*KKS1, 1, 1, 2*KS1, 1, K2KS3, 0, SN3) 
                  SN3 = SN3*DSQRT(DBLE((K2KS3 + 1)*(2*KKS1 + 1))) 
                  IF (MOD(KKS1,2) /= 0) SN3 = -SN3 
                  CAS1 = CAS1 + SN1*SN3*RECS1*RECS2*RECS3 
               ENDIF 
               IF (ISES2*ISES4 == 0) CYCLE  
               CALL SIXJ (2*KKS1, 1, 1, 2*KS2, 1, K2KS3, 0, SN4) 
               SN4 = SN4*DSQRT(DBLE((K2KS3 + 1)*(2*KKS1 + 1))) 
               IF (MOD(KKS1,2) /= 0) SN4 = -SN4 
               CBS1 = CBS1 + SN2*SN4*RECS1*RECS2*RECS3 
            END DO 
         END DO 
         IF (IREZ == 1) THEN 
            IF (MOD(KS2,2) /= 0) THEN 
               CAS1 = -CAS1 
               CBS1 = -CBS1 
            ENDIF 
         ELSE 
            IF (MOD(KS1,2) /= 0) THEN 
               CAS1 = -CAS1 
               CBS1 = -CBS1 
            ENDIF 
         ENDIF 
         RW1(KS1+1,KS2+8) = CAS1 
         RW1(KS1+1,KS2+10) = CBS1 
      ELSE 
         CAS1 = RW1(KS1+1,KS2+8) 
         CBS1 = RW1(KS1+1,KS2+10) 
      ENDIF 
      ISES1 = 1 
      ISES2 = 1 
      IF (DABS(CAS1) < EPS) ISES1 = 0 
      IF (DABS(CBS1) < EPS) ISES2 = 0 
      IF (ISES1 + ISES2 == 0) RETURN  
      IF (IWAA(2,1,KL1+1,KL2+1) == 0) THEN 
         IWAA(2,1,KL1+1,KL2+1) = 1 
         DO I3 = IP1, IG1 
            KKL3 = I3 - 1 
            IP2 = ITREXG(LJ(IA),LJ(IB),KKL3,LJ(IC),IKK2) + 1 
            IF (IKK2 == 0) CYCLE  
            IG2 = IP2 + IKK2 - 1 
            ISES1L = IXJTIK(2*KL2,2*KL1,2*K,2*KKL3,2*ID4(3),2*ID2(3)) 
            ISES2L = IXJTIK(2*KL1,2*KL2,2*K,2*KKL3,2*ID4(3),2*ID2(3)) 
            IF (ISES1L + ISES2L == 0) CYCLE  
            IF (IW2(1,I3) == 0) THEN 
               IW2(1,I3) = 1 
               CALL RLSP4B (2, IC, ID, 2*KKL3, 2*ID4(3), 2*K, 0, IAT, R) 
               IF (IAT == 0) THEN 
                  RECL1 = ZERO 
               ELSE 
                  CALL RLSP4B (2, IC, ID, 2*KKL3, 2*ID4(3), 2*K, 1, IAT, RECL1) 
               ENDIF 
               RW2(1,I3) = RECL1 
            ELSE 
               RECL1 = RW2(1,I3) 
            ENDIF 
            IF (DABS(RECL1) < EPS) CYCLE  
            IF (ISES1L /= 0) CALL SIXJ (2*KL2, 2*KL1, 2*K, 2*KKL3, 2*ID4(3), 2*&
               ID2(3), 0, S3) 
            IF (ISES2L /= 0) CALL SIXJ (2*KL1, 2*KL2, 2*K, 2*KKL3, 2*ID4(3), 2*&
               ID2(3), 0, S4) 
            DO I4 = IP2, IG2 
               KKL1 = I4 - 1 
               ISES3L = IXJTIK(2*KL1,2*ID3(3),2*ID1(3),2*KKL1,2*ID2(3),2*KKL3) 
               ISES4L = IXJTIK(2*KL2,2*ID3(3),2*ID1(3),2*KKL1,2*ID2(3),2*KKL3) 
               IF (ISES3L + ISES4L == 0) CYCLE  
               IF (IRL(I4,I3) == 0) THEN 
                  IRL(I4,I3) = 1 
                  IAT = 0 
                  CALL DLSA4 (2, IB, IC, 2*KKL1, 2*ID3(3), 2*KKL3, 0, IAT, R) 
                  IF (IAT == 0) THEN 
                     RECL3 = ZERO 
                  ELSE 
                     CALL DLSA4 (2, IB, IC, 2*KKL1, 2*ID3(3), 2*KKL3, 1, IAT, &
                        RECL3) 
                  ENDIF 
                  RL(I4,I3) = RECL3 
               ELSE 
                  RECL3 = RL(I4,I3) 
               ENDIF 
               IF (DABS(RECL3) <= EPS) CYCLE  
               IF (IW2(2,I4) == 0) THEN 
                  IW2(2,I4) = 1 
                  CALL RLSP4A (2, IA, IB, IC, 2*ID1(3), 2*ID2(3), 2*KKL1, 0, &
                     IAT, R) 
                  IF (IAT == 0) THEN 
                     RECL2 = ZERO 
                  ELSE 
                     CALL RLSP4A (2, IA, IB, IC, 2*ID1(3), 2*ID2(3), 2*KKL1, 1&
                        , IAT, RECL2) 
                  ENDIF 
                  RW2(2,I4) = RECL2 
               ELSE 
                  RECL2 = RW2(2,I4) 
               ENDIF 
               IF (DABS(RECL2) <= EPS) CYCLE  
               IF (ISES1L*ISES3L /= 0) THEN 
                  CALL SIXJ (2*KL1, 2*ID3(3), 2*ID1(3), 2*KKL1, 2*ID2(3), 2*&
                     KKL3, 0, SLA1) 
                  SLA1 = SLA1*DSQRT(DBLE((2*KKL1 + 1)*(2*KKL3 + 1))) 
                  IF (MOD(KKL1,2) /= 0) SLA1 = -SLA1 
                  CA1 = CA1 + S3*SLA1*RECL1*RECL2*RECL3 
               ENDIF 
               IF (ISES2L*ISES4L == 0) CYCLE  
               CALL SIXJ (2*KL2, 2*ID3(3), 2*ID1(3), 2*KKL1, 2*ID2(3), 2*KKL3, &
                  0, SLB1) 
               SLB1 = SLB1*DSQRT(DBLE((2*KKL1 + 1)*(2*KKL3 + 1))) 
               IF (MOD(KKL1,2) /= 0) SLB1 = -SLB1 
               CB1 = CB1 + S4*SLB1*RECL1*RECL2*RECL3 
            END DO 
         END DO 
         IF (IREZ == 1) THEN 
            IF (MOD(ID2(3)+ID3(3)+KL2,2) /= 0) THEN 
               CA1 = -CA1 
               CB1 = -CB1 
            ENDIF 
         ELSE 
            IF (MOD(ID1(3)+ID4(3)+KL1,2) /= 0) THEN 
               CA1 = -CA1 
               CB1 = -CB1 
            ENDIF 
         ENDIF 
         RWAA(2,1,KL1+1,KL2+1) = CA1 
         RWAA(2,2,KL1+1,KL2+1) = CB1 
      ELSE 
         CA1 = RWAA(2,1,KL1+1,KL2+1) 
         CB1 = RWAA(2,2,KL1+1,KL2+1) 
      ENDIF 
      CA1 = -CA1*CAS1*HALF*RW1(1,1) 
      CB1 = -CB1*CBS1*HALF*RW1(1,1) 
      RETURN  
      END SUBROUTINE TWO56 
