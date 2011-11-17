!
!     -------------------------------------------------------------
!      N O N R E L A T 3 1
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
!                                              N'2 = N2 + 1        *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!
      SUBROUTINE NONRELAT31(IA, IB, IIA, IIB, IIC, IID, IIRE) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE MEDEFN_C 
      USE DIAGNL_C 
      USE OPERAT_C 
      USE TRK_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:17:07  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE recoupls0_I 
      USE calculation_I 
      USE recoupls2_I 
      USE hibff_I 
      USE a1a2ls_I 
      USE savenon_I 
      USE itrexg_I 
      USE coulombls_I 
      USE orbitorbit_I 
      USE ixjtik_I 
      USE a1aw2ls_I 
      USE sixj_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: IA 
      INTEGER  :: IB 
      INTEGER , INTENT(IN) :: IIA 
      INTEGER , INTENT(IN) :: IIB 
      INTEGER , INTENT(IN) :: IIC 
      INTEGER , INTENT(IN) :: IID 
      INTEGER , INTENT(IN) :: IIRE 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IAA, IBB, IAT, LIA2, LIB2, LA, LB, NN, IB1, II, LC, LD, IP2, &
         IKK, IG2, I2, KL, I3, L12, I13, JS12, IFAZ, KLL 
      REAL(DOUBLE) :: REC, RECS, RECL, QM1, QM2, A1, A2, WW, B, AB, AC, AD, SI&
         , AA, ABB 
!-----------------------------------------------
      IF (IHSH <= 1) RETURN  
      IF (IA == IB) RETURN  
      IF (IA < IB) THEN 
         IAA = IA 
         IBB = IB 
      ELSE 
         IAA = IB 
         IBB = IA 
      ENDIF 
      IF (.NOT.RECOUPLS0(1,IAA,IBB,IBB,IBB,0)) RETURN  
      IF (.NOT.RECOUPLS0(2,IAA,IBB,IBB,IBB,1)) RETURN  
      IF (.NOT.RECOUPLS0(3,IAA,IBB,IBB,IBB,1)) RETURN  
      CALL RECOUPLS2 (3, IAA, IBB, 1, 0, IAT, REC) 
      IF (IAT == 0) RETURN  
      LIA2 = LJ(IA)*2 
      LIB2 = LJ(IB)*2 
      CALL RECOUPLS2 (2, IAA, IBB, LIB2, 0, IAT, REC) 
      IF (IAT == 0) RETURN  
      CALL RECOUPLS2 (3, IAA, IBB, 1, 1, IAT, RECS) 
      CALL RECOUPLS2 (2, IAA, IBB, LIB2, 1, IAT, RECL) 
      QM1 = HALF 
      QM2 = -HALF 
      A1 = ZERO 
      A2 = ZERO 
      LA = IJFUL(IAA) 
      LB = IJFUL(IBB) 
      CALL HIBFF (IA, IB, IA, IA, 2) 
      IF (LJ(IA) == LJ(IB)) THEN 
         CALL A1A2LS (IK1, IK2, BK1, BK2, ID1, ID2, BD1, BD2, QM2, QM1, WW) 
         IF (DABS(WW) > EPS) THEN 
            B = WW*RECL*RECS*HALF*DSQRT(DBLE(4*LJ(IA)+2)) 
            NN = 0 
            IB1 = IBB - 1 
            NN = SUM(NOSH1(IAA:IB1)) 
            IF ((NN/2)*2 == NN) B = -B 
            IF (DABS(B) > EPS) CALL SAVENON (4, B, 0, 0, LA, 0, LB, JA, JB, 0) 
         ENDIF 
      ENDIF 
      IF (IIRE == 0) RETURN  
!
!     CASES 2111   + + - -        TRANSFORM TO  1112   + - - +
!           1211                                1112
!
      LA = IJFUL(IIA) 
      LB = IJFUL(IIB) 
      LC = IJFUL(IIC) 
      LD = IJFUL(IID) 
      IP2 = ITREXG(LJ(IA),LJ(IA),LJ(IA),LJ(IB),IKK) + 1 
      IF (IKK <= 0) RETURN  
      IG2 = IP2 + IKK - 1 
      DO I2 = IP2, IG2 
         KL = I2 - 1 
         IF (.NOT.CALCULATION(KL)) CYCLE  
         IF (ICOLOM + ISOTOP == 1) CALL COULOMBLS (LJ(IB), LJ(IA), LJ(IA), LJ(&
            IA), KL, A1) 
         IF (IORBORB == 1) CALL ORBITORBIT (LJ(IB), LJ(IA), LJ(IA), LJ(IA), KL&
             + 1, A2) 
         IF (DABS(A1) + DABS(A2) <= EPS) CYCLE  
         AB = ZERO 
         DO I3 = IP2, IG2 
            L12 = I3 - 1 
            IF (IXJTIK(LIB2,LIA2,KL*2,LIA2,LIA2,L12*2) == 0) CYCLE  
            AC = ZERO 
            DO I13 = 1, 2 
               JS12 = I13 - 1 
               IFAZ = L12 + JS12 
               IF ((IFAZ/2)*2 /= IFAZ) CYCLE  
               CALL A1AW2LS (IK2, IK1, BK2, BK1, ID2, ID1, BD2, BD1, L12, JS12&
                  , QM1, QM1, QM2, QM2, AD) 
               AD = AD*DSQRT(DBLE(2*JS12 + 1)) 
               AC = AC + AD 
            END DO 
            CALL SIXJ (LIB2, LIA2, KL*2, LIA2, LIA2, L12*2, 0, SI) 
            AA = AC*SI*DSQRT(DBLE(2*L12 + 1)) 
            IFAZ = KL + L12 
            IF ((IFAZ/2)*2 /= IFAZ) AA = -AA 
            AB = AB + AA 
         END DO 
         AB = -AB*RECL*RECS 
         IF (DABS(AB) <= EPS) CYCLE  
         NN = 0 
         IB1 = IBB - 1 
         NN = SUM(NOSH1(IAA:IB1)) 
         IF ((NN/2)*2 == NN) AB = -AB 
         ABB = AB 
         IF (ICOLOM + ISOTOP == 1) THEN 
            IF (DABS(A1) > EPS) THEN 
               AB = A1*AB 
               CALL SAVENON (3, AB, KL, LA, LB, LC, LD, JA, JB, 0) 
            ENDIF 
         ENDIF 
! Orbit 2111
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
      END SUBROUTINE NONRELAT31 
