!
!     -------------------------------------------------------------
!      N O N R E L A T 3 3
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 2313, 3231, 3213, 2331    *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 - 1   *
!                                                   N'2 = N2 + 1   *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!
      SUBROUTINE NONRELAT33(IA, IB, IC, IREZ, IIA, IIB, IIC, IID) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE MEDEFN_C 
      USE DIAGNL_C 
      USE OPERAT_C 
      USE TRK_C 
      USE TRK2_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:17:07  11/16/01  
!...Switches:                     
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE recoupls0_I 
      USE calculation_I 
      USE eile_I 
      USE hibff_I 
      USE itrexg_I 
      USE recoupls3_I 
      USE a1a2w3ls_I 
      USE jfaze_I 
      USE coulombls_I 
      USE savenon_I 
      USE orbitorbit_I 
      USE ixjtik_I 
      USE sixj_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: IA 
      INTEGER  :: IB 
      INTEGER  :: IC 
      INTEGER , INTENT(IN) :: IREZ 
      INTEGER , INTENT(IN) :: IIA 
      INTEGER , INTENT(IN) :: IIB 
      INTEGER , INTENT(IN) :: IIC 
      INTEGER , INTENT(IN) :: IID 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: LMAX = 10 
      INTEGER, PARAMETER :: L2MAX = 2*LMAX 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(2) :: IATS 
      INTEGER :: IAA, IBB, ICC, LA, LB, LC, LD, IP1, IKK, IG1, IAT, LIA2, LIB2&
         , LIC2, I4, KL, IFAZP, NN, IB1, II, I1, KL2, II1, IP2, IG2, I2, I3, &
         L12, I13, KLL 
      REAL(DOUBLE), DIMENSION(L2MAX,2) :: PM 
      REAL(DOUBLE) :: QM1, QM2, A1, A2, RECS1, RECS2, REC, RAG1, RAG2, AA, AB, &
         AC, AD, SI, ABB 
!-----------------------------------------------
!
      IF (IHSH <= 2) RETURN  
      CALL EILE (IA, IB, IC, IAA, IBB, ICC) 
      IF (.NOT.RECOUPLS0(1,IAA,ICC,IBB,IBB,0)) RETURN  
      IF (.NOT.RECOUPLS0(2,IAA,ICC,IBB,IBB,2)) RETURN  
      IF (.NOT.RECOUPLS0(3,IAA,ICC,IBB,IBB,2)) RETURN  
      LA = IJFUL(IIA) 
      LB = IJFUL(IIB) 
      LC = IJFUL(IIC) 
      LD = IJFUL(IID) 
      QM1 = HALF 
      QM2 = -HALF 
      A1 = ZERO 
      A2 = ZERO 
      CALL HIBFF (IA, IB, IC, IA, 3) 
      IP1 = ITREXG(LJ(IB),LJ(IA),LJ(IC),LJ(IC),IKK) + 1 
      IF (IKK <= 0) RETURN  
      IG1 = IP1 + IKK - 1 
      RECS1 = ZERO 
      RECS2 = ZERO 
      CALL RECOUPLS3 (3, IB, IA, IC, 1, 1, 0, 0, IATS(1), REC) 
      IF (IATS(1) /= 0) CALL RECOUPLS3 (3, IB, IA, IC, 1, 1, 0, 1, IAT, RECS1) 
      CALL RECOUPLS3 (3, IB, IA, IC, 1, 1, 2, 0, IATS(2), REC) 
      IF (IATS(2) /= 0) CALL RECOUPLS3 (3, IB, IA, IC, 1, 1, 2, 1, IAT, RECS2) 
      IF (IATS(1) + IATS(2) == 0) RETURN  
      LIA2 = LJ(IA)*2 
      LIB2 = LJ(IB)*2 
      LIC2 = LJ(IC)*2 
      DO I4 = IP1, IG1 
         KL = I4 - 1 
         PM(I4,1) = ZERO 
         PM(I4,2) = ZERO 
         CALL RECOUPLS3 (2, IB, IA, IC, LIB2, LIA2, KL*2, 0, IAT, REC) 
         IF (IAT == 0) CYCLE  
         CALL A1A2W3LS (IK2, IK1, IK3, BK2, BK1, BK3, ID2, ID1, ID3, BD2, BD1, &
            BD3, KL, 0, QM1, QM2, QM1, QM2, RAG1) 
         CALL A1A2W3LS (IK2, IK1, IK3, BK2, BK1, BK3, ID2, ID1, ID3, BD2, BD1, &
            BD3, KL, 1, QM1, QM2, QM1, QM2, RAG2) 
         IF (DABS(RAG1) + DABS(RAG2) <= EPS) CYCLE  
         CALL RECOUPLS3 (2, IB, IA, IC, LIB2, LIA2, KL*2, 1, IAT, REC) 
         PM(I4,1) = REC*RECS1*RAG1 
         PM(I4,2) = REC*RECS2*RAG2 
      END DO 
      IFAZP = JFAZE(IB,IA,IC,IC) 
      IF (IA < IB) THEN 
         IAA = IA 
         IBB = IB 
      ELSE 
         IAA = IB 
         IBB = IA 
      ENDIF 
      NN = 0 
      IB1 = IBB - 1 
      NN = SUM(NOSH1(IAA:IB1)) 
      IF ((NN/2)*2 == NN) IFAZP = -IFAZP 
      IF (IREZ == 2) GO TO 5 
!
!     CASES 2313   + + - -        TRANSFORM TO  2133   + - + -
!           3231                                2133
!
    6 CONTINUE 
      IF (IATS(1) /= 0) THEN 
         DO I1 = IP1, IG1 
            KL = I1 - 1 
            IF (.NOT.CALCULATION(KL)) CYCLE  
            IF (ICOLOM + ISOTOP == 1) THEN 
               CALL COULOMBLS (LJ(IB), LJ(IC), LJ(IA), LJ(IC), KL, AA) 
               AA = AA*PM(I1,1)/DSQRT(DBLE(2*KL + 1)) 
               AA = AA*TWO*DBLE(IFAZP) 
               IF (DABS(AA) > EPS) CALL SAVENON (3, AA, KL, LA, LB, LC, LD, JA&
                  , JB, 0) 
            ENDIF 
! Orbit 2313
            IF (IORBORB /= 1) CYCLE  
            KL2 = KL + 2 
            CALL ORBITORBIT (LJ(IB), LJ(IC), LJ(IA), LJ(IC), KL2, AA) 
            IF (DABS(AA) <= EPS) CYCLE  
            II1 = I1 + 1 
            AA = AA*PM(II1,1)/DBLE(2*(KL2 - 1) + 1) 
            AA = AA*DBLE(IFAZP) 
            IF (DABS(AA) > EPS) THEN 
               CALL SAVENON (9, AA, KL, LA, LB, LC, LD, JA, JB, 0) 
               CALL SAVENON (9, AA, KL, LB, LA, LD, LC, JA, JB, 0) 
            ENDIF 
!                IF(DABS(AA).GT.EPS)
!     :                         WRITE(79,556) AA,KL,LJ(IA),LJ(IB),JA,JB
  556       FORMAT(1X,'1221 O-O','AA=',F17.7,'K=',I3,'LA=',2I3,'JA=',I4,'JB=',I&
               4) 
         END DO 
      ENDIF 
      IF (IREZ == 2) GO TO 7 
!
!     CASES 3213   + + - -        TRANSFORM TO  2133   + - + -
!           2331                                2133
!
    5 CONTINUE 
      IP2 = ITREXG(LJ(IC),LJ(IA),LJ(IB),LJ(IC),IKK) + 1 
      IF (IKK > 0) THEN 
         IG2 = IP2 + IKK - 1 
         DO I2 = IP2, IG2 
            KL = I2 - 1 
            IF (.NOT.CALCULATION(KL)) CYCLE  
            IF (ICOLOM + ISOTOP == 1) CALL COULOMBLS (LJ(IC), LJ(IB), LJ(IA), &
               LJ(IC), KL, A1) 
            IF (IORBORB == 1) CALL ORBITORBIT (LJ(IC), LJ(IB), LJ(IA), LJ(IC), &
               KL + 1, A2) 
            IF (DABS(A1) + DABS(A2) <= EPS) CYCLE  
            AB = ZERO 
            DO I3 = IP1, IG1 
               L12 = I3 - 1 
               IF (IXJTIK(LIA2,LIC2,KL*2,LIC2,LIB2,L12*2) == 0) CYCLE  
               AC = ZERO 
               DO I13 = 1, 2 
                  IF (IATS(I13) == 0) CYCLE  
                  AD = PM(I3,I13)*DSQRT(DBLE(2*(I13 - 1) + 1)) 
                  IF (I13 == 1) AD = -AD 
                  AC = AC + AD 
               END DO 
               IF (DABS(AC) <= EPS) CYCLE  
               CALL SIXJ (LIA2, LIC2, KL*2, LIC2, LIB2, L12*2, 0, SI) 
               AA = AC*SI*SQRT(DBLE(2*L12 + 1)) 
               AB = AB + AA 
            END DO 
            IF (DABS(AB) <= EPS) CYCLE  
            AB = AB*DBLE(IFAZP) 
            ABB = AB 
            IF (ICOLOM + ISOTOP == 1) THEN 
               IF (DABS(A1) > EPS) THEN 
                  AB = A1*AB 
                  CALL SAVENON (3, AB, KL, LA, LB, LD, LC, JA, JB, 0) 
               ENDIF 
            ENDIF 
! Orbit 3213
            IF (IORBORB /= 1) CYCLE  
            IF (DABS(A2) <= EPS) CYCLE  
            KLL = KL - 1 
            ABB = A2*ABB*HALF/DSQRT(DBLE(2*KL + 1)) 
            CALL SAVENON (9, ABB, KLL, LA, LB, LD, LC, JA, JB, 0) 
            CALL SAVENON (9, ABB, KLL, LB, LA, LC, LD, JA, JB, 0) 
!                 WRITE(79,656) ABB,KLL,LJ(IA),LJ(IB),JA,JB
  656       FORMAT(1X,'1221 O-O','AA=',F17.7,'K=',I3,'LA=',2I3,'JA=',I4,'JB=',I&
               4) 
         END DO 
      ENDIF 
      IF (IREZ == 2) GO TO 6 
    7 CONTINUE 
      RETURN  
      END SUBROUTINE NONRELAT33 
