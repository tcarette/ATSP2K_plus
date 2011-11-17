!
!     -------------------------------------------------------------
!      N O N R E L A T 1
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :           N'1 = N1        *
!                                                  N'2 = N2        *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!                                                                  *
!     Universite Libre de Bruxelles, Brussels                      *
!                                                  December 1995   *
!
      SUBROUTINE NONRELAT1(IA, IB, IIRE) 
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
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE recoupls0_I 
      USE calculation_I 
      USE hibff_I 
      USE w1_I 
      USE savenon_I 
      USE coulombls_I 
      USE wwls1_I 
      USE avera_I 
      USE orbitorbit_I 
      USE ittk_I 
      USE recoupls2_I 
      USE w1w2ls_I 
      USE ixjtik_I 
      USE sixj_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: IA 
      INTEGER  :: IB 
      INTEGER , INTENT(IN) :: IIRE 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: LMAX = 10 
      INTEGER, PARAMETER :: L2MAX = 2*LMAX 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LIA2, LA, IN, IG2, I2, KL, KL2, KL1, LIB2, LB, IG1, I4, IAT, &
         I1, II1, IP2, IGAL, I3, L12, I13, KLL 
      REAL(DOUBLE), DIMENSION(L2MAX,2) :: PMG, RAG 
      REAL(DOUBLE) :: QM1, QM2, A1, A2, A, RECOPL, B, AA, RE, REP, RA, AB, AC, &
         SI, ABB 
      LOGICAL , DIMENSION(2) :: IATT 
      LOGICAL :: IATTT 
!-----------------------------------------------
!
      LIA2 = LJ(IA)*2 
      LA = IJFUL(IA) 
      QM1 = HALF 
      QM2 = -HALF 
      A1 = ZERO 
      A2 = ZERO 
      IF (IA == IB) THEN 
!
!     THE CASE 1111   + + - -
!
         IF (JA /= JB) THEN 
            DO IN = 1, 3 
               IF (RECOUPLS0(IN,IA,IA,IA,IA,0)) CYCLE  
               RETURN  
            END DO 
         ENDIF 
         CALL HIBFF (IA, IA, IA, IA, 1) 
         CALL W1 (IK1, BK1, ID1, BD1, 0, 0, QM1, QM2, A) 
         RECOPL = ONE/DSQRT(DBLE(IK1(6)+1)*(IK1(5)+1)) 
         IF (ABS(A) > EPS) THEN 
            B = A*RECOPL*HALF*DSQRT(DBLE(4*LJ(IA)+2)) 
!BAY            IF(DABS(B).GT.EPS) CALL SAVENON(4,B,0,0,LA,0,LA,JA,JB,0)
            IF (DABS(B) > EPS) THEN 
               IF (JA /= JB) THEN 
                  WRITE (*, '(A,I5,A,I5)') 'Configurations', JA, 'and', JB, &
                     'are identical' 
                  STOP  
               ELSE 
                  CALL SAVENON (4, B, 0, 0, LA, 0, LA, JA, JB, 0) 
               ENDIF 
            ENDIF 
            A = A/DSQRT(DBLE(4*LJ(IA)+2)) 
         ENDIF 
         IF (IIRE == 0) RETURN  
         IF (ICOLOM + IORBORB /= 0) THEN 
            IG2 = LIA2 + 1 
            DO I2 = 1, IG2, 2 
               KL = I2 - 1 
               IF (ICOLOM == 1) THEN 
!GGf                                with    f-shell   ***  beginning
                  IF (JA /= JB) THEN 
!GGf                                with    f-shell   ***  e n d
                     CALL COULOMBLS (LJ(IA), LJ(IA), LJ(IA), LJ(IA), KL, A1) 
                     IF (DABS(A1) > EPS) THEN 
                        CALL WWLS1 (IK1, BK1, ID1, BD1, KL, 0, QM1, QM2, QM1, &
                           QM2, AA) 
                        AA = AA/DSQRT(DBLE(2*KL + 1)) + A 
                        AA = AA*A1*RECOPL 
                        IF (DABS(AA) > EPS) CALL SAVENON (1, AA, KL, 0, LA, 0, &
                           LA, JA, JB, 0) 
                     ENDIF 
                  ELSE IF (LJ(IA) > 3) THEN 
                     CALL COULOMBLS (LJ(IA), LJ(IA), LJ(IA), LJ(IA), KL, A1) 
                     IF (DABS(A1) > EPS) THEN 
                        CALL WWLS1 (IK1, BK1, ID1, BD1, KL, 0, QM1, QM2, QM1, &
                           QM2, AA) 
                        AA = AA/DSQRT(DBLE(2*KL + 1)) + A 
                        AA = AA*A1*RECOPL 
                        IF (DABS(AA) > EPS) CALL SAVENON (1, AA, KL, 0, LA, 0, &
                           LA, JA, JB, 0) 
                     ENDIF 
!GGf                                with    f-shell   ***  beginning
                  ELSE 
                     CALL AVERA (ID1(3), ID1(1), KL, ID1(4), AA) 
                     IF (DABS(AA) > EPS) CALL SAVENON (1, AA, KL, 0, LA, 0, LA&
                        , JA, JB, 0) 
!GGf                                with    f-shell   ***  e n d
                  ENDIF 
               ENDIF 
! Orbit 1111
               IF (IORBORB /= 1) CYCLE  
               KL2 = KL + 2 
               KL1 = KL + 1 
               CALL ORBITORBIT (LJ(IA), LJ(IA), LJ(IA), LJ(IA), KL2, A1) 
               IF (DABS(A1) <= EPS) CYCLE  
               CALL WWLS1 (IK1, BK1, ID1, BD1, KL1, 0, QM1, QM2, QM1, QM2, AA) 
               AA = AA/DBLE(2*KL1 + 1) - A/DSQRT(DBLE(2*KL1 + 1)) 
               AA = AA*A1*RECOPL 
               IF (DABS(AA) > EPS) CALL SAVENON (9, AA, KL, LA, LA, LA, LA, JA&
                  , JB, 0) 
!                IF(DABS(AA).GT.EPS)
!     :                               WRITE(77,555) AA,KL,LA,JA,JB
  555          FORMAT(1X,'1111 O-O','AA=',F17.7,'K=',I3,'LA=',I3,'J=',2I4) 
            END DO 
            RETURN  
         ENDIF 
      ENDIF 
      IF (IIRE == 0) RETURN  
      IF (IHSH <= 1) RETURN  
      IF (ISOTOP == 1) THEN 
         IF (ITTK(LJ(IA),LJ(IB),1) == 0) RETURN  
      ENDIF 
      IATT(1) = .TRUE. 
      IATT(2) = .TRUE. 
      IF (JA /= JB) THEN 
         IF (.NOT.RECOUPLS0(1,IA,IB,IB,IB,0)) RETURN  
         IF (.NOT.RECOUPLS0(2,IA,IB,IB,IB,1)) RETURN  
         IATT(1) = RECOUPLS0(3,IA,IB,IB,IB,0) 
         IATT(2) = RECOUPLS0(3,IA,IB,IB,IB,1) 
         IATTT = .TRUE. 
         IF (IATT(1) .OR. IATT(2)) IATTT = .FALSE. 
         IF (IATTT) RETURN  
      ENDIF 
      CALL HIBFF (IA, IB, IA, IA, 2) 
      LIB2 = LJ(IB)*2 
      LB = IJFUL(IB) 
      IG1 = MIN(LIA2,LIB2) + 1 
      DO I4 = 1, IG1 
         KL = I4 - 1 
         RAG(I4,1) = ZERO 
         PMG(I4,1) = ZERO 
         RAG(I4,2) = ZERO 
         PMG(I4,2) = ZERO 
         CALL RECOUPLS2 (2, IA, IB, KL*2, 0, IAT, RE) 
         IF (IAT == 0) CYCLE  
         CALL RECOUPLS2 (2, IA, IB, KL*2, 1, IAT, RE) 
         IF (IATT(1)) THEN 
            CALL RECOUPLS2 (3, IA, IB, 0, 0, IAT, REP) 
            IF (IAT /= 0) THEN 
               CALL W1W2LS (KL, 0, KL, 0, QM1, QM2, QM1, QM2, RA) 
               RAG(I4,1) = RA 
               PMG(I4,1) = RE/DSQRT(DBLE((IK1(6)+1)*(IK2(6)+1))) 
            ENDIF 
         ENDIF 
         IF (.NOT.IATT(2)) CYCLE  
         CALL RECOUPLS2 (3, IA, IB, 2, 0, IAT, REP) 
         IF (IAT == 0) CYCLE  
         CALL W1W2LS (KL, 1, KL, 1, QM1, QM2, QM1, QM2, RA) 
         IF (DABS(RA) <= EPS) CYCLE  
         RAG(I4,2) = RA 
         CALL RECOUPLS2 (3, IA, IB, 2, 1, IAT, REP) 
         PMG(I4,2) = RE*REP 
      END DO 
!
!     CASES 1212   + + - -        TRANSFORM TO  1122   + - + -
!           2121                                1122
!
      IF (IATT(1)) THEN 
         IF (ICOLOM + IORBORB /= 0) THEN 
            DO I1 = 1, IG1, 2 
               KL = I1 - 1 
               IF (ICOLOM == 1) THEN 
                  CALL COULOMBLS (LJ(IA), LJ(IB), LJ(IA), LJ(IB), KL, AA) 
                  IF (DABS(AA) >= EPS) THEN 
                     AA = AA*PMG(I1,1) 
                     AA = AA*RAG(I1,1) 
                     AA = AA*TWO/DSQRT(DBLE(2*KL + 1)) 
                     IF (DABS(AA) > EPS) CALL SAVENON (1, AA, KL, 0, LA, 0, LB&
                        , JA, JB, 0) 
                  ENDIF 
               ENDIF 
! Orbit 1212
               IF (IORBORB /= 1) CYCLE  
               KL2 = KL + 2 
               CALL ORBITORBIT (LJ(IA), LJ(IB), LJ(IA), LJ(IB), KL2, AA) 
               IF (DABS(AA) <= EPS) CYCLE  
               II1 = I1 + 1 
               AA = AA*PMG(II1,1) 
               AA = AA*RAG(II1,1) 
               AA = AA/DBLE(2*(KL2 - 1) + 1) 
               IF (DABS(AA) > EPS) THEN 
                  CALL SAVENON (9, AA, KL, LA, LB, LA, LB, JA, JB, 0) 
                  CALL SAVENON (9, AA, KL, LB, LA, LB, LA, JA, JB, 0) 
               ENDIF 
!                IF(DABS(AA).GT.EPS)
!     :                         WRITE(78,556) AA,KL,LJ(IA),LJ(IB),JA,JB
  556          FORMAT(1X,'1212 O-O','AA=',F17.7,'K=',I3,'LA=',2I3,'JA=',I4,&
                  'JB=',I4) 
            END DO 
         ENDIF 
      ENDIF 
!
!     CASES 1221   + + - -        TRANSFORM TO  1122   + - + -
!           2112                                1122
!
      IP2 = IABS(LJ(IB)-LJ(IA)) + 1 
      IG2 = LJ(IB) + LJ(IA) + 1 
      IGAL = 2 
      IF (IORBORB == 1) IGAL = 1 
      DO I2 = IP2, IG2, IGAL 
         KL = I2 - 1 
         IF (.NOT.CALCULATION(KL)) CYCLE  
         IF (ICOLOM + ISOTOP == 1) CALL COULOMBLS (LJ(IA), LJ(IB), LJ(IB), LJ(&
            IA), KL, A1) 
         IF (IORBORB == 1) CALL ORBITORBIT (LJ(IA), LJ(IB), LJ(IB), LJ(IA), KL&
             + 1, A2) 
         IF (DABS(A1) + DABS(A2) <= EPS) CYCLE  
         AB = ZERO 
         DO I3 = 1, IG1 
            L12 = I3 - 1 
            AC = ZERO 
            IF (IXJTIK(LIA2,LIB2,KL*2,LIB2,LIA2,L12*2) == 0) CYCLE  
            DO I13 = 1, 2 
               IF (.NOT.IATT(I13)) CYCLE  
               AA = PMG(I3,I13) 
               AA = AA*RAG(I3,I13) 
               AA = AA*DSQRT(DBLE(2*(I13 - 1) + 1)) 
               IF (I13 == 1) AA = -AA 
               AC = AC + AA 
            END DO 
            IF (DABS(AC) <= EPS) CYCLE  
            CALL SIXJ (LIA2, LIB2, KL*2, LIB2, LIA2, L12*2, 0, SI) 
            AA = AC*SI*DSQRT(DBLE(2*L12 + 1)) 
            AB = AB + AA 
         END DO 
         IF (DABS(AB) <= EPS) CYCLE  
         ABB = AB 
         IF (ICOLOM + ISOTOP == 1) THEN 
            AB = A1*AB 
            IF (DABS(AB) > EPS) CALL SAVENON (2, AB, KL, 0, LA, 0, LB, JA, JB, &
               0) 
         ENDIF 
! Orbit 1221
         IF (IORBORB /= 1) CYCLE  
         IF (DABS(A2) <= EPS) CYCLE  
         KLL = KL - 1 
         ABB = A2*ABB*HALF/DSQRT(DBLE(2*KL + 1)) 
         CALL SAVENON (9, ABB, KLL, LA, LB, LB, LA, JA, JB, 0) 
         CALL SAVENON (9, ABB, KLL, LB, LA, LA, LB, JA, JB, 0) 
!                 WRITE(79,656) ABB,KLL,LJ(IA),LJ(IB),JA,JB
  656    FORMAT(1X,'1221 O-O','AA=',F17.7,'K=',I3,'LA=',2I3,'JA=',I4,'JB=',I4) 
      END DO 
      RETURN  
      END SUBROUTINE NONRELAT1 
