!
!     -------------------------------------------------------------
!      N O N R E L A T 5 2
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 1324, 3142, 1342, 3124    *
!                                                   ( IREZ = 1),   *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 - 1   *
!                                                   N'2 = N2 + 1   *
!                                                   N'3 = N3 - 1   *
!                                                   N'4 = N4 + 1   *
!     AND    2413, 4231, 2431, 4213                 ( IREZ = 2),   *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 + 1   *
!                                                   N'2 = N2 - 1   *
!                                                   N'3 = N3 + 1   *
!                                                   N'4 = N4 - 1   *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!
      SUBROUTINE NONRELAT52(IA, IB, IC, ID, IREZ) 
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
      USE hibff_I 
      USE a1a2a3a4ls_I 
      USE recoupls4_I 
      USE itrexg_I 
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
      INTEGER  :: ID 
      INTEGER , INTENT(IN) :: IREZ 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: LMAX = 10 
      INTEGER, PARAMETER :: L2MAX = 2*LMAX 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(2) :: IATS 
      INTEGER , DIMENSION(L2MAX) :: IATL 
      INTEGER :: I13, JS12, IAT, IP1, IKK, IG1, LA, LB, LC, LD, LIA2, LIB2, &
         LIC2, LID2, IAT1, I3, L12, IFAZP, NN, IB1, II, ID11, I2, KL, IFAZ, KL2&
         , II2, IP2, IG2, KLL 
      REAL(DOUBLE), DIMENSION(2) :: PMS 
      REAL(DOUBLE), DIMENSION(L2MAX) :: PML 
      REAL(DOUBLE) :: QM1, QM2, A1, A2, ANG, REC, RECS, RECL, AB, AC, AA, SI, &
         ABB 
!-----------------------------------------------
!
      IF (IHSH <= 3) RETURN  
      IF (.NOT.RECOUPLS0(1,IA,ID,IC,IB,0)) RETURN  
      IF (.NOT.RECOUPLS0(2,IA,ID,IC,IB,3)) RETURN  
      IF (.NOT.RECOUPLS0(3,IA,ID,IC,IB,3)) RETURN  
      IF (IREZ == 1) THEN 
         QM1 = HALF 
         QM2 = -HALF 
      ELSE 
         QM1 = -HALF 
         QM2 = HALF 
      ENDIF 
      A1 = ZERO 
      A2 = ZERO 
      CALL HIBFF (IA, IB, IC, ID, 4) 
      CALL A1A2A3A4LS (IK1, IK2, IK3, IK4, BK1, BK2, BK3, BK4, ID1, ID2, ID3, &
         ID4, BD1, BD2, BD3, BD4, QM1, QM2, QM1, QM2, ANG) 
      IF (DABS(ANG) < EPS) RETURN  
      PMS(1) = ZERO 
      PMS(2) = ZERO 
      DO I13 = 1, 2 
         JS12 = I13 - 1 
         CALL RECOUPLS4 (3, IA, IB, IC, ID, 1, 1, 1, 1, JS12*2, 0, IAT, REC) 
         IATS(I13) = IAT 
         IF (IATS(I13) == 0) CYCLE  
         CALL RECOUPLS4 (3, IA, IB, IC, ID, 1, 1, 1, 1, JS12*2, 1, IAT, RECS) 
         PMS(I13) = RECS 
      END DO 
      IF (IATS(1) + IATS(2) == 0) RETURN  
      IP1 = ITREXG(LJ(IB),LJ(IA),LJ(IC),LJ(ID),IKK) + 1 
      IF (IKK <= 0) RETURN  
      IG1 = IP1 + IKK - 1 
      LA = IJFUL(IA) 
      LB = IJFUL(IB) 
      LC = IJFUL(IC) 
      LD = IJFUL(ID) 
      LIA2 = LJ(IA)*2 
      LIB2 = LJ(IB)*2 
      LIC2 = LJ(IC)*2 
      LID2 = LJ(ID)*2 
      IAT1 = 0 
      DO I3 = IP1, IG1 
         PML(I3) = ZERO 
         L12 = I3 - 1 
         CALL RECOUPLS4 (2, IA, IB, IC, ID, LIA2, LIB2, LIC2, LID2, L12*2, 0, &
            IAT, REC) 
         IATL(I3) = IAT 
         IF (IATL(I3) /= 0) THEN 
            CALL RECOUPLS4 (2, IA, IB, IC, ID, LIA2, LIB2, LIC2, LID2, L12*2, 1&
               , IAT, RECL) 
            PML(I3) = RECL 
         ENDIF 
         IAT1 = IAT1 + IATL(I3) 
      END DO 
      IF (IAT1 == 0) RETURN  
      IFAZP = 1 
      NN = 0 
      IB1 = IB - 1 
      NN = SUM(NOSH1(IA:IB1)) 
      IF ((NN/2)*2 == NN) IFAZP = -IFAZP 
      NN = 0 
      ID11 = ID - 1 
      NN = SUM(NOSH1(IC:ID11)) 
      IF ((NN/2)*2 == NN) IFAZP = -IFAZP 
!
!     CASES 1324   + + - -        TRANSFORM TO  1234   + - + -
!           1342                                1234
!                                                    (IREZ = 1)
!     OR
!     CASES 2413   + + - -        TRANSFORM TO  1234   - + - +
!           4231                                1234
!                                                    (IREZ = 2)
      IF (IATS(1) /= 0) THEN 
         DO I2 = IP1, IG1 
            KL = I2 - 1 
            IF (.NOT.CALCULATION(KL)) CYCLE  
            IF (ICOLOM + ISOTOP == 1) THEN 
               IF (IATL(I2) /= 0) THEN 
                  CALL COULOMBLS (LJ(IA), LJ(IC), LJ(IB), LJ(ID), KL, A1) 
                  IF (DABS(A1) > EPS) THEN 
                     AB = TWO*A1*ANG*PML(I2)*PMS(1)*DBLE(IFAZP)/DSQRT(DBLE(2*KL&
                         + 1)) 
                     IFAZ = 2 
                     IF (IREZ == 2) IFAZ = IK1(3) + IK4(3) + IK2(3) + IK3(3) 
                     IF ((IFAZ/2)*2 /= IFAZ) AB = -AB 
                     IF (DABS(AB) > EPS) THEN 
                        IF (IREZ == 1) CALL SAVENON (3, AB, KL, LA, LC, LB, LD&
                           , JA, JB, 0) 
                        IF (IREZ == 2) CALL SAVENON (3, AB, KL, LB, LD, LA, LC&
                           , JA, JB, 0) 
                     ENDIF 
                  ENDIF 
               ENDIF 
            ENDIF 
! Orbit 1324
            IF (IORBORB /= 1) CYCLE  
            KL2 = KL + 2 
            CALL ORBITORBIT (LJ(IA), LJ(IC), LJ(IB), LJ(ID), KL2, A2) 
            IF (DABS(A2) <= EPS) CYCLE  
            II2 = I2 + 1 
            AB = A2*ANG*PML(II2)*PMS(1)*DBLE(IFAZP)/DBLE(2*(KL2 - 1) + 1) 
            IF (DABS(AB) > EPS) THEN 
               IFAZ = 2 
               IF (IREZ == 2) IFAZ = IK1(3) + IK4(3) + IK2(3) + IK3(3) 
               IF ((IFAZ/2)*2 /= IFAZ) AB = -AB 
               CALL SAVENON (9, AB, KL, LA, LC, LB, LD, JA, JB, 0) 
               CALL SAVENON (9, AB, KL, LC, LA, LD, LB, JA, JB, 0) 
            ENDIF 
!                IF(DABS(AB).GT.EPS)
!     :                         WRITE(79,556) AB,KL,LJ(IA),LJ(IB),JA,JB
  556       FORMAT(1X,'1221 O-O','AA=',F17.7,'K=',I3,'LA=',2I3,'JA=',I4,'JB=',I&
               4) 
         END DO 
      ENDIF 
!
!     CASES 1342   + + - -        TRANSFORM TO  1234   + - + -
!           3124                                1234
!                                                    (IREZ = 1)
!     OR
!     CASES 2431   + + - -        TRANSFORM TO  1234   - + - +
!           4213                                1234
!                                                    (IREZ = 2)
      IP2 = ITREXG(LJ(IA),LJ(ID),LJ(IB),LJ(IC),IKK) + 1 
      IF (IKK <= 0) RETURN  
      IG2 = IP2 + IKK - 1 
      DO I2 = IP2, IG2 
         KL = I2 - 1 
         IF (.NOT.CALCULATION(KL)) CYCLE  
         IF (ICOLOM + ISOTOP == 1) CALL COULOMBLS (LJ(IA), LJ(IC), LJ(ID), LJ(&
            IB), KL, A1) 
         IF (IORBORB == 1) CALL ORBITORBIT (LJ(IA), LJ(IC), LJ(ID), LJ(IB), KL&
             + 1, A2) 
         IF (DABS(A1) + DABS(A2) <= EPS) CYCLE  
         AB = ZERO 
         DO I3 = IP1, IG1 
            L12 = I3 - 1 
            IF (IATL(I3) == 0) CYCLE  
            IF (IXJTIK(LIA2,LID2,KL*2,LIC2,LIB2,L12*2) == 0) CYCLE  
            AC = ZERO 
            DO I13 = 1, 2 
               IF (IATS(I13) == 0) CYCLE  
               JS12 = I13 - 1 
               AA = PMS(I13)*DSQRT(DBLE(2*JS12 + 1)) 
               IF (I13 == 1) AA = -AA 
               AC = AC + AA 
            END DO 
            CALL SIXJ (LIA2, LID2, KL*2, LIC2, LIB2, L12*2, 0, SI) 
            AA = AC*PML(I3)*SI*DSQRT(DBLE(2*L12 + 1)) 
            AB = AB + AA 
         END DO 
         IF (DABS(AB) <= EPS) CYCLE  
         AB = ANG*AB*DBLE(IFAZP) 
         IFAZ = 2 
         IF (IREZ == 2) IFAZ = IK1(3) - IK4(3) + IK2(3) + IK3(3) 
         IF ((IFAZ/2)*2 /= IFAZ) AB = -AB 
         ABB = AB 
         IF (ICOLOM + ISOTOP == 1) THEN 
            IF (DABS(A1) > EPS) THEN 
               AB = A1*AB 
               IF (IREZ == 1) CALL SAVENON (3, AB, KL, LA, LC, LD, LB, JA, JB, &
                  0) 
               IF (IREZ == 2) CALL SAVENON (3, AB, KL, LB, LD, LC, LA, JA, JB, &
                  0) 
            ENDIF 
         ENDIF 
! Orbit 1342
         IF (IORBORB /= 1) CYCLE  
         IF (DABS(A2) <= EPS) CYCLE  
         KLL = KL - 1 
         ABB = A2*ABB*HALF/DSQRT(DBLE(2*KL + 1)) 
         CALL SAVENON (9, ABB, KLL, LA, LC, LD, LB, JA, JB, 0) 
         CALL SAVENON (9, ABB, KLL, LC, LA, LB, LD, JA, JB, 0) 
!                 WRITE(79,655) ABB,KLL,LJ(IA),LJ(IB),JA,JB
  655    FORMAT(1X,'1221 O-O','AA=',F17.7,'K=',I3,'LA=',2I3,'JA=',I4,'JB=',I4) 
      END DO 
      RETURN  
      END SUBROUTINE NONRELAT52 
