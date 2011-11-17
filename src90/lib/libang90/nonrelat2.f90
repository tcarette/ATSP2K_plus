!
!     -------------------------------------------------------------
!      N O N R E L A T 2
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 2        *
!                                              N'2 = N2 + 2        *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!
      SUBROUTINE NONRELAT2(IA, IB) 
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
      USE ittk_I 
      USE hibff_I 
      USE recoupls2_I 
      USE coulombls_I 
      USE orbitorbit_I 
      USE ixjtik_I 
      USE w1w2ls_I 
      USE sixj_I 
      USE savenon_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: IA 
      INTEGER  :: IB 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IAA, IBB, I13, JS12, IAT, LA, LB, LIA2, LIB2, IG1, IP2, IG2, &
         IGAL, I2, KL, I3, L12, IFAZ, KLL 
      REAL(DOUBLE), DIMENSION(2) :: PMS 
      REAL(DOUBLE) :: QM1, QM2, A1, A2, REC, RECS, AB, AC, AA, RECL, SI, ABB 
      LOGICAL , DIMENSION(2) :: IATT 
      LOGICAL :: IATTT 
!-----------------------------------------------
      IF (IHSH <= 1) RETURN  
      IF (IA == IB) RETURN  
      IF (ISOTOP == 1) THEN 
         IF (ITTK(LJ(IA),LJ(IB),1) == 0) RETURN  
      ENDIF 
      IF (IA < IB) THEN 
         IAA = IA 
         IBB = IB 
      ELSE 
         IAA = IB 
         IBB = IA 
      ENDIF 
      IF (.NOT.RECOUPLS0(1,IAA,IBB,IBB,IBB,0)) RETURN  
      IF (.NOT.RECOUPLS0(2,IAA,IBB,IBB,IBB,1)) RETURN  
      IATT(1) = RECOUPLS0(3,IAA,IBB,IBB,IBB,0) 
      IATT(2) = RECOUPLS0(3,IAA,IBB,IBB,IBB,1) 
      IATTT = .TRUE. 
      IF (IATT(1) .OR. IATT(2)) IATTT = .FALSE. 
      IF (IATTT) RETURN  
      QM1 = HALF 
      QM2 = -HALF 
      A1 = ZERO 
      A2 = ZERO 
      CALL HIBFF (IA, IB, IA, IA, 2) 
!
!     THE CASE 1122   + + - -
!
      PMS(1) = ZERO 
      PMS(2) = ZERO 
      DO I13 = 1, 2 
         IF (.NOT.IATT(I13)) CYCLE  
         JS12 = I13 - 1 
         CALL RECOUPLS2 (3, IAA, IBB, JS12*2, 0, IAT, REC) 
         IF (IAT == 0) CYCLE  
         CALL RECOUPLS2 (3, IAA, IBB, JS12*2, 1, IAT, RECS) 
         PMS(I13) = RECS 
      END DO 
      LA = IJFUL(IA) 
      LB = IJFUL(IB) 
      LIA2 = LJ(IA)*2 
      LIB2 = LJ(IB)*2 
      IG1 = MIN(LIA2,LIB2) + 1 
      IP2 = IABS(LJ(IB)-LJ(IA)) + 1 
      IG2 = LJ(IB) + LJ(IA) + 1 
      IGAL = 2 
      IF (IORBORB == 1) IGAL = 1 
      DO I2 = IP2, IG2, IGAL 
         KL = I2 - 1 
         IF (.NOT.CALCULATION(KL)) CYCLE  
         IF (ICOLOM + ISOTOP == 1) CALL COULOMBLS (LJ(IA), LJ(IA), LJ(IB), LJ(&
            IB), KL, A1) 
         IF (IORBORB == 1) CALL ORBITORBIT (LJ(IA), LJ(IA), LJ(IB), LJ(IB), KL&
             + 1, A2) 
         IF (DABS(A1) + DABS(A2) <= EPS) CYCLE  
         A1 = -HALF*A1 
         AB = ZERO 
         DO I3 = 1, IG1 
            L12 = I3 - 1 
            CALL RECOUPLS2 (2, IAA, IBB, L12*2, 0, IAT, REC) 
            IF (IAT == 0) CYCLE  
            IF (IXJTIK(LIA2,LIB2,KL*2,LIB2,LIA2,L12*2) == 0) CYCLE  
            AC = ZERO 
            DO I13 = 1, 2 
               IF (.NOT.IATT(I13)) CYCLE  
               JS12 = I13 - 1 
               CALL W1W2LS (L12, JS12, L12, JS12, QM1, QM1, QM2, QM2, AA) 
               AA = AA*PMS(I13)*DSQRT(DBLE(2*JS12 + 1)) 
               AC = AC + AA 
            END DO 
            CALL RECOUPLS2 (2, IAA, IBB, L12*2, 1, IAT, RECL) 
            CALL SIXJ (LIA2, LIB2, KL*2, LIB2, LIA2, L12*2, 0, SI) 
            AA = AC*RECL*SI*DSQRT(DBLE(2*L12 + 1)) 
            IFAZ = IK1(3) + IK2(3) + KL + L12 
            IF ((IFAZ/2)*2 /= IFAZ) AA = -AA 
            AB = AB + AA 
         END DO 
         IF (DABS(AB) <= EPS) CYCLE  
         ABB = AB 
         IF (ICOLOM + ISOTOP == 1) THEN 
            AB = A1*AB 
            IF (DABS(AB) > EPS) CALL SAVENON (3, AB, KL, LA, LA, LB, LB, JA, JB&
               , 0) 
         ENDIF 
! Orbit 1122
         IF (IORBORB /= 1) CYCLE  
         IF (DABS(A2) <= EPS) CYCLE  
         KLL = KL - 1 
         ABB = -A2*ABB*HALF/DSQRT(DBLE(2*KL + 1)) 
         CALL SAVENON (9, ABB, KLL, LA, LA, LB, LB, JA, JB, 0) 
!                 WRITE(79,656) ABB,KLL,LJ(IA),LJ(IB),JA,JB
  656    FORMAT(1X,'1221 O-O','AA=',F17.7,'K=',I3,'LA=',2I3,'JA=',I4,'JB=',I4) 
      END DO 
      RETURN  
      END SUBROUTINE NONRELAT2 
