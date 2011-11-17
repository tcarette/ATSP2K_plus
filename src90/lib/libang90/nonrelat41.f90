!
!     -------------------------------------------------------------
!      N O N R E L A T 4 1
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 + 1        *
!                                              N'2 = N2 + 1        *
!                                              N'3 = N3 - 2,       *
!     WHEN IREZ = 1   ........................
!                                              N'1 = N1 - 1        *
!                                              N'2 = N2 - 1        *
!                                              N'3 = N3 + 2,       *
!     WHEN IREZ = 2   ........................
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!
      SUBROUTINE NONRELAT41(IA, IB, IC, IREZ, IIA, IIB, IIC, IID) 
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
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE recoupls0_I 
      USE calculation_I 
      USE eile_I 
      USE hibff_I 
      USE recoupls3_I 
      USE itrexg_I 
      USE jfaze_I 
      USE coulombls_I 
      USE orbitorbit_I 
      USE ixjtik_I 
      USE a1a2w3ls_I 
      USE sixj_I 
      USE savenon_I 
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
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(2) :: IATS 
      INTEGER :: IAA, IBB, ICC, I13, JS12, IAT, LA, LB, LC, LD, LIA2, LIB2, &
         LIC2, IP1, IKK, IG1, IP2, IG2, IFAZP, NN, IB1, II, I2, KL, I3, L12, &
         IFAZ, KLL 
      REAL(DOUBLE), DIMENSION(2) :: PMS 
      REAL(DOUBLE) :: QM1, QM2, A1, A2, REC, RECS, AB, AC, AA, RECL, SI, ABB 
!-----------------------------------------------
      IF (IHSH <= 2) RETURN  
      CALL EILE (IA, IB, IC, IAA, IBB, ICC) 
      IF (.NOT.RECOUPLS0(1,IAA,ICC,IBB,IBB,0)) RETURN  
      IF (.NOT.RECOUPLS0(2,IAA,ICC,IBB,IBB,2)) RETURN  
      IF (.NOT.RECOUPLS0(3,IAA,ICC,IBB,IBB,2)) RETURN  
      IF (IREZ == 1) THEN 
         QM1 = -HALF 
         QM2 = HALF 
      ELSE 
         QM1 = HALF 
         QM2 = -HALF 
      ENDIF 
      A1 = ZERO 
      A2 = ZERO 
      CALL HIBFF (IA, IB, IC, IA, 3) 
!
!     CASES 3312   + + - -        TRANSFORM TO  1233   - - + +
!           3321                                1233
!                                                    (IREZ = 1)
!     OR
!     CASES 1233   + + - -        TRANSFORM TO  1233   + + - -
!           2133                                1233
!                                                    (IREZ = 2)
      PMS(1) = ZERO 
      PMS(2) = ZERO 
      DO I13 = 1, 2 
         JS12 = I13 - 1 
         CALL RECOUPLS3 (3, IA, IB, IC, 1, 1, JS12*2, 0, IAT, REC) 
         IATS(I13) = IAT 
         IF (IATS(I13) == 0) CYCLE  
         CALL RECOUPLS3 (3, IA, IB, IC, 1, 1, JS12*2, 1, IAT, RECS) 
         PMS(I13) = RECS 
      END DO 
      IF (IATS(1) + IATS(2) == 0) RETURN  
      LA = IJFUL(IIA) 
      LB = IJFUL(IIB) 
      LC = IJFUL(IIC) 
      LD = IJFUL(IID) 
      LIA2 = LJ(IA)*2 
      LIB2 = LJ(IB)*2 
      LIC2 = LJ(IC)*2 
      IP1 = ITREXG(LJ(IB),LJ(IA),LJ(IC),LJ(IC),IKK) + 1 
      IF (IKK <= 0) RETURN  
      IG1 = IP1 + IKK - 1 
      IP2 = ITREXG(LJ(IC),LJ(IA),LJ(IB),LJ(IC),IKK) + 1 
      IF (IKK <= 0) RETURN  
      IG2 = IP2 + IKK - 1 
      IFAZP = JFAZE(IC,IA,IB,IC) 
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
      DO I2 = IP2, IG2 
         KL = I2 - 1 
         IF (.NOT.CALCULATION(KL)) CYCLE  
         IF (ICOLOM + ISOTOP == 1) CALL COULOMBLS (LJ(IA), LJ(IB), LJ(IC), LJ(&
            IC), KL, A1) 
         IF (IORBORB == 1) CALL ORBITORBIT (LJ(IA), LJ(IB), LJ(IC), LJ(IC), KL&
             + 1, A2) 
         IF (DABS(A1) + DABS(A2) <= EPS) CYCLE  
         AB = ZERO 
         DO I3 = IP1, IG1 
            L12 = I3 - 1 
            CALL RECOUPLS3 (2, IA, IB, IC, LIA2, LIB2, L12*2, 0, IAT, REC) 
            IF (IAT == 0) CYCLE  
            IF (IXJTIK(LIC2,LIB2,KL*2,LIA2,LIC2,L12*2) == 0) CYCLE  
            AC = ZERO 
            DO I13 = 1, 2 
               IF (IATS(I13) == 0) CYCLE  
               JS12 = I13 - 1 
               IFAZ = L12 + JS12 
               IF ((IFAZ/2)*2 /= IFAZ) CYCLE  
               CALL A1A2W3LS (IK1, IK2, IK3, BK1, BK2, BK3, ID1, ID2, ID3, BD1&
                  , BD2, BD3, L12, JS12, QM1, QM1, QM2, QM2, AA) 
               AA = AA*PMS(I13)*DSQRT(DBLE(2*JS12 + 1)) 
               AC = AC + AA 
            END DO 
            CALL RECOUPLS3 (2, IA, IB, IC, LIA2, LIB2, L12*2, 1, IAT, RECL) 
            CALL SIXJ (LIC2, LIB2, KL*2, LIA2, LIC2, L12*2, 0, SI) 
            AA = AC*RECL*SI*DSQRT(DBLE(2*L12 + 1)) 
            IFAZ = IK1(3) + IK3(3) + KL + L12 
            IF (IREZ == 2) IFAZ = IK2(3) + IK3(3) + KL + L12 
            IF ((IFAZ/2)*2 /= IFAZ) AA = -AA 
            AB = AB + AA 
         END DO 
         IF (DABS(AB) <= EPS) CYCLE  
         AB = -AB*DBLE(IFAZP) 
         ABB = AB 
         IF (ICOLOM + ISOTOP == 1) THEN 
            IF (DABS(A1) > EPS) THEN 
               AB = A1*AB 
               CALL SAVENON (3, AB, KL, LA, LB, LC, LD, JA, JB, 0) 
            ENDIF 
         ENDIF 
! Orbit 3312
         IF (IORBORB /= 1) CYCLE  
         IF (DABS(A2) <= EPS) CYCLE  
         KLL = KL - 1 
         ABB = A2*ABB*HALF/DSQRT(DBLE(2*KL + 1)) 
         CALL SAVENON (9, ABB, KLL, LA, LB, LC, LD, JA, JB, 0) 
         CALL SAVENON (9, ABB, KLL, LB, LA, LD, LC, JA, JB, 0) 
!                 WRITE(79,656) ABB,KLL,LJ(IA),LJ(IB),JA,JB
  656    FORMAT(1X,'1221 O-O','AA=',F17.7,'K=',I3,'LA=',2I3,'JA=',I4,'JB=',I4) 
      END DO 
      RETURN  
      END SUBROUTINE NONRELAT41 
