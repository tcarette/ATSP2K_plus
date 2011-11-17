!
!     -------------------------------------------------------------
!      T W O 5 5
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
!                                              N'2 = N2 + 1        *
!     ( IREZ = 1 )                             N'3 = N3 - 1        *
!                                              N'4 = N4 + 1        *
!                                                                  *
!     CASES 1423   + + - -        TRANSFORM TO  1234     + - - +   *
!           4132                                1234               *
!                                                                  *
!                                                                  *
!                                                                  *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 + 1        *
!                                              N'2 = N2 - 1        *
!     ( IREZ = 2 )                             N'3 = N3 + 1        *
!                                              N'4 = N4 - 1        *
!                                                                  *
!     CASES 2314   + + - -        TRANSFORM TO  1234     - + + -   *
!           3241                                1234               *
!                                                                  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE TWO55(KL1, KS1, KL2, KS2, K, IA, IB, IC, ID, IREZ, CA1, CB1) 
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
      USE rlsp4a_I 
      USE ixjtik_I 
      USE rlsp4b_I 
      USE dlsa4_I 
      USE sixj_I 
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
      INTEGER :: JN, JI, JIS, IPS1, IKK1, IGS1, IPS2, IP1, IG1, IAT, ISES3, &
         ISES4, ISES1, ISES2, K2KS3, KIIS1, ISES1L, ISES2L, I3, KKL3 
      REAL(DOUBLE) :: CAS1, CBS1, R, RECS2A, RECS2B, RECS1, RECS3A, RECS3B, SN1&
         , SN2, RECL2A, RECL2B, RECL1, RECL3A, RECL3B, S3, S4 
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
      IF (IW1(KS1+1,KS2+4) == 0) THEN 
         IW1(KS1+1,KS2+4) = 1 
         CAS1 = ZERO 
         CBS1 = ZERO 
         IF (IW1(KS1+1,3) == 0) THEN 
            IW1(KS1+1,3) = 1 
            CALL RLSP4A (3, IA, IB, IC, 1, 1, 2*KS1, 0, IAT, R) 
            IF (IAT == 0) THEN 
               RECS2A = ZERO 
            ELSE 
               CALL RLSP4A (3, IA, IB, IC, 1, 1, 2*KS1, 1, IAT, RECS2A) 
            ENDIF 
            RW1(KS1+1,3) = RECS2A 
         ELSE 
            RECS2A = RW1(KS1+1,3) 
         ENDIF 
         IF (IW1(KS2+1,3) == 0) THEN 
            IW1(KS2+1,3) = 1 
            CALL RLSP4A (3, IA, IB, IC, 1, 1, 2*KS2, 0, IAT, R) 
            IF (IAT == 0) THEN 
               RECS2B = ZERO 
            ELSE 
               CALL RLSP4A (3, IA, IB, IC, 1, 1, 2*KS2, 1, IAT, RECS2B) 
            ENDIF 
            RW1(KS2+1,3) = RECS2B 
         ELSE 
            RECS2B = RW1(KS2+1,3) 
         ENDIF 
         ISES3 = 1 
         ISES4 = 1 
         IF (DABS(RECS2A) < EPS) ISES3 = 0 
         IF (DABS(RECS2B) < EPS) ISES4 = 0 
         IF (ISES3 + ISES4 /= 0) THEN 
            ISES1 = 0 
            ISES2 = 0 
            DO K2KS3 = IPS1, IGS1, 2 
               IF (ISES3 /= 0) ISES1 = IXJTIK(2*KS2,2*KS1,2*K,K2KS3,1,1) 
               IF (ISES4 /= 0) ISES2 = IXJTIK(2*KS1,2*KS2,2*K,K2KS3,1,1) 
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
               IF (DABS(RECS1) > EPS) THEN 
                  IF (ISES1 /= 0) THEN 
                     IF (IRS(KS1+1,KIIS1) == 0) THEN 
                        IRS(KS1+1,KIIS1) = 1 
                        IAT = 0 
                        CALL DLSA4 (3, IB, IC, 2*KS1, 1, K2KS3, 0, IAT, R) 
                        IF (IAT == 0) THEN 
                           RECS3A = ZERO 
                        ELSE 
                           CALL DLSA4 (3, IB, IC, 2*KS1, 1, K2KS3, 1, IAT, &
                              RECS3A) 
                        ENDIF 
                        RS(KS1+1,KIIS1) = RECS3A 
                     ELSE 
                        RECS3A = RS(KS1+1,KIIS1) 
                     ENDIF 
                     IF (DABS(RECS3A) < EPS) ISES1 = 0 
                  ENDIF 
                  IF (ISES2 /= 0) THEN 
                     IF (IRS(KS2+1,KIIS1) == 0) THEN 
                        IRS(KS2+1,KIIS1) = 1 
                        IAT = 0 
                        CALL DLSA4 (3, IB, IC, 2*KS2, 1, K2KS3, 0, IAT, R) 
                        IF (IAT == 0) THEN 
                           RECS3B = ZERO 
                        ELSE 
                           CALL DLSA4 (3, IB, IC, 2*KS2, 1, K2KS3, 1, IAT, &
                              RECS3B) 
                        ENDIF 
                        RS(KS2+1,KIIS1) = RECS3B 
                     ELSE 
                        RECS3B = RS(KS2+1,KIIS1) 
                     ENDIF 
                     IF (DABS(RECS3B) < EPS) ISES2 = 0 
                  ENDIF 
               ENDIF 
               IF (ISES1 + ISES2 == 0) CYCLE  
               IF (ISES1 /= 0) THEN 
                  CALL SIXJ (2*KS2, 2*KS1, 2*K, K2KS3, 1, 1, 0, SN1) 
                  SN1 = SN1*DSQRT(DBLE(K2KS3 + 1)) 
                  CAS1 = CAS1 + SN1*RECS1*RECS2A*RECS3A 
               ENDIF 
               IF (ISES2 == 0) CYCLE  
               CALL SIXJ (2*KS1, 2*KS2, 2*K, K2KS3, 1, 1, 0, SN2) 
               SN2 = SN2*DSQRT(DBLE(K2KS3 + 1)) 
               CBS1 = CBS1 + SN2*RECS1*RECS2B*RECS3B 
            END DO 
         ENDIF 
         CAS1 = CAS1/DSQRT(DBLE(2*KS1 + 1)) 
         CBS1 = CBS1/DSQRT(DBLE(2*KS2 + 1)) 
         IF (IREZ == 1) THEN 
            IF (MOD(KS1 + KS2,2) /= 0) CAS1 = -CAS1 
         ELSE 
            IF (MOD(KS1 + KS2,2) /= 0) CBS1 = -CBS1 
         ENDIF 
         RW1(KS1+1,KS2+4) = CAS1 
         RW1(KS1+1,KS2+6) = CBS1 
      ELSE 
         CAS1 = RW1(KS1+1,KS2+4) 
         CBS1 = RW1(KS1+1,KS2+6) 
      ENDIF 
      ISES1 = 1 
      ISES2 = 1 
      IF (DABS(CAS1) < EPS) ISES1 = 0 
      IF (DABS(CBS1) < EPS) ISES2 = 0 
      IF (ISES1 + ISES2 == 0) RETURN  
      IF (IWAA(1,1,KL1+1,KL2+1) == 0) THEN 
         IWAA(1,1,KL1+1,KL2+1) = 1 
         IF (IW2(2,KL1+1) == 0) THEN 
            IW2(2,KL1+1) = 1 
            CALL RLSP4A (2, IA, IB, IC, 2*ID1(3), 2*ID2(3), 2*KL1, 0, IAT, R) 
            IF (IAT == 0) THEN 
               RECL2A = ZERO 
            ELSE 
               CALL RLSP4A (2, IA, IB, IC, 2*ID1(3), 2*ID2(3), 2*KL1, 1, IAT, &
                  RECL2A) 
            ENDIF 
            RW2(2,KL1+1) = RECL2A 
         ELSE 
            RECL2A = RW2(2,KL1+1) 
         ENDIF 
         IF (IW2(2,KL2+1) == 0) THEN 
            IW2(2,KL2+1) = 1 
            CALL RLSP4A (2, IA, IB, IC, 2*ID1(3), 2*ID2(3), 2*KL2, 0, IAT, R) 
            IF (IAT == 0) THEN 
               RECL2B = ZERO 
            ELSE 
               CALL RLSP4A (2, IA, IB, IC, 2*ID1(3), 2*ID2(3), 2*KL2, 1, IAT, &
                  RECL2B) 
            ENDIF 
            RW2(2,KL2+1) = RECL2B 
         ELSE 
            RECL2B = RW2(2,KL2+1) 
         ENDIF 
         ISES3 = 1 
         ISES4 = 1 
         IF (DABS(RECL2A) < EPS) ISES3 = 0 
         IF (DABS(RECL2B) < EPS) ISES4 = 0 
         IF (ISES3 + ISES4 /= 0) THEN 
            ISES1L = 0 
            ISES2L = 0 
            DO I3 = IP1, IG1 
               KKL3 = I3 - 1 
               IF (ISES3 /= 0) ISES1L = IXJTIK(2*KL2,2*KL1,2*K,2*KKL3,2*ID4(3),&
                  2*ID3(3)) 
               IF (ISES4 /= 0) ISES2L = IXJTIK(2*KL1,2*KL2,2*K,2*KKL3,2*ID4(3),&
                  2*ID3(3)) 
               IF (ISES1L + ISES2L == 0) CYCLE  
               IF (IW2(1,I3) == 0) THEN 
                  IW2(1,I3) = 1 
                  CALL RLSP4B (2, IC, ID, 2*KKL3, 2*ID4(3), 2*K, 0, IAT, R) 
                  IF (IAT == 0) THEN 
                     RECL1 = ZERO 
                  ELSE 
                     CALL RLSP4B (2, IC, ID, 2*KKL3, 2*ID4(3), 2*K, 1, IAT, &
                        RECL1) 
                  ENDIF 
                  RW2(1,I3) = RECL1 
               ELSE 
                  RECL1 = RW2(1,I3) 
               ENDIF 
               IF (DABS(RECL1) <= EPS) CYCLE  
               IF (ISES1L /= 0) THEN 
                  IF (IRL(KL1+1,I3) == 0) THEN 
                     IRL(KL1+1,I3) = 1 
                     IAT = 0 
                     CALL DLSA4 (2, IB, IC, 2*KL1, 2*ID3(3), 2*KKL3, 0, IAT, R) 
                     IF (IAT == 0) THEN 
                        RECL3A = ZERO 
                     ELSE 
                        CALL DLSA4 (2, IB, IC, 2*KL1, 2*ID3(3), 2*KKL3, 1, IAT&
                           , RECL3A) 
                     ENDIF 
                     RL(KL1+1,I3) = RECL3A 
                  ELSE 
                     RECL3A = RL(KL1+1,I3) 
                  ENDIF 
                  IF (DABS(RECL3A) < EPS) ISES1L = 0 
               ELSE 
                  ISES1L = 0 
               ENDIF 
               IF (ISES2L /= 0) THEN 
                  IF (IRL(KL2+1,I3) == 0) THEN 
                     IRL(KL2+1,I3) = 1 
                     IAT = 0 
                     CALL DLSA4 (2, IB, IC, 2*KL2, 2*ID3(3), 2*KKL3, 0, IAT, R) 
                     IF (IAT == 0) THEN 
                        RECL3B = ZERO 
                     ELSE 
                        CALL DLSA4 (2, IB, IC, 2*KL2, 2*ID3(3), 2*KKL3, 1, IAT&
                           , RECL3B) 
                     ENDIF 
                     RL(KL2+1,I3) = RECL3B 
                  ELSE 
                     RECL3B = RL(KL2+1,I3) 
                  ENDIF 
                  IF (DABS(RECL3B) < EPS) ISES2L = 0 
               ELSE 
                  ISES2L = 0 
               ENDIF 
               IF (ISES1L + ISES2L == 0) CYCLE  
               IF (ISES1L /= 0) THEN 
                  CALL SIXJ (2*KL2, 2*KL1, 2*K, 2*KKL3, 2*ID4(3), 2*ID3(3), 0, &
                     S3) 
                  CA1 = CA1 + S3*RECL1*RECL2A*RECL3A*DSQRT(DBLE(2*KKL3 + 1)) 
               ENDIF 
               IF (ISES2L == 0) CYCLE  
               CALL SIXJ (2*KL1, 2*KL2, 2*K, 2*KKL3, 2*ID4(3), 2*ID3(3), 0, S4) 
               CB1 = CB1 + S4*RECL1*RECL2B*RECL3B*DSQRT(DBLE(2*KKL3 + 1)) 
            END DO 
         ENDIF 
         CA1 = CA1/DSQRT(DBLE(2*KL1 + 1)) 
         CB1 = CB1/DSQRT(DBLE(2*KL2 + 1)) 
         IF (IREZ == 1) THEN 
            IF (MOD(KL1 + KL2 + 1,2) /= 0) CA1 = -CA1 
            CB1 = -CB1 
         ELSE 
            IF (MOD(ID1(3)+ID2(3)+ID3(3)+ID4(3)+1,2) /= 0) CA1 = -CA1 
            IF (MOD(ID1(3)+ID2(3)+ID3(3)+ID4(3)+KL1+KL2+1,2) /= 0) CB1 = -CB1 
         ENDIF 
         RWAA(1,1,KL1+1,KL2+1) = CA1 
         RWAA(1,2,KL1+1,KL2+1) = CB1 
      ELSE 
         CA1 = RWAA(1,1,KL1+1,KL2+1) 
         CB1 = RWAA(1,2,KL1+1,KL2+1) 
      ENDIF 
      CA1 = HALF*CA1*CAS1*RW1(1,1) 
      CB1 = HALF*CB1*CBS1*RW1(1,1) 
      RETURN  
      END SUBROUTINE TWO55 
