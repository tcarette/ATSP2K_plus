!                                                                  :
!     ..........................................................   :
!        ii)   two - particle operator                             :
!     ..........................................................   :
!
!     -------------------------------------------------------------
!      T W O 1
!     -------------------------------------------------------------
!                                                                  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
!                                                                  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE TWO1(IIA, IIB, IIRE) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE MEDEFN_C 
      USE TRK_C 
      USE PERMAT_C 
      USE KAMPAS_C 
      USE CASEOP_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:23:13  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE two13_I 
      USE ittk_I 
      USE hibff_I 
      USE rlsp0_I 
      USE rlsp1_I 
      USE ss1111_I 
      USE soo1111_I 
      USE ss1212_I 
      USE soo1212_I 
      USE ss1221_I 
      USE sooc_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IIA 
      INTEGER , INTENT(IN) :: IIB 
      INTEGER , INTENT(IN) :: IIRE 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KGAL, IA, IAT, LIA2, IG, I, KL, IB, I1, I2, LIB2, IP 
      REAL(DOUBLE) :: REC, RECS, RECL 
!-----------------------------------------------
!
      IF (IOCASE == 1) THEN 
! Spin - spin      operator
         KGAL = 2 
      ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
         KGAL = 1 
      ENDIF 
!
      IF (IIA == IIB) THEN 
         IA = IIA 
         IF (ITTK(J1QN1(IA,2)-1,J1QN2(IA,2)-1,2*KGAL) == 0) RETURN  
         IF (ITTK(J1QN1(IA,3)-1,J1QN2(IA,3)-1,2*KGAL) == 0) RETURN  
         CALL HIBFF (IA, IA, IA, IA, 1) 
         IF (ID1(4) < 2) RETURN  
         CALL RLSP0 (1, IA, IA, 0, IAT) 
         IF (IAT == 0) RETURN  
         CALL RLSP0 (2, IA, IA, 2*KGAL, IAT) 
         IF (IAT == 0) RETURN  
         CALL RLSP0 (3, IA, IA, 2*KGAL, IAT) 
         IF (IAT == 0) RETURN  
         CALL RLSP1 (3, IA, 2*KGAL, 0, IAT, REC) 
         IF (IAT == 0) RETURN  
         CALL RLSP1 (2, IA, 2*KGAL, 0, IAT, REC) 
         IF (IAT == 0) RETURN  
         CALL RLSP1 (3, IA, 2*KGAL, 1, IAT, RECS) 
         IF (DABS(RECS) < EPS) RETURN  
         CALL RLSP1 (2, IA, 2*KGAL, 1, IAT, RECL) 
         IF (DABS(RECL) < EPS) RETURN  
         RS(1,1) = RECS 
         RL(1,1) = RECL 
         LIA2 = LJ(IA)*2 
         IG = LIA2 + 1 
         IRS(1,1) = 0 
         DO I = 1, IG, 2 
            KL = I - 1 
            IF (IOCASE == 1) THEN 
! Spin - spin      operator
               CALL SS1111 (IG, KL, IA) 
            ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
               CALL SOO1111 (IG, KL, IA) 
!            CALL SOO1111P(IG,KL,IA)
            ENDIF 
         END DO 
      ELSE 
         IF (IIRE == 0) RETURN  
         IF (IHSH <= 1) RETURN  
         IA = MIN0(IIA,IIB) 
         IB = MAX0(IIA,IIB) 
         CALL RLSP0 (1, IA, IB, 0, IAT) 
         IF (IAT == 0) RETURN  
         CALL RLSP0 (2, IA, IB, 2*KGAL, IAT) 
         IF (IAT == 0) RETURN  
         CALL RLSP0 (3, IA, IB, 2*KGAL, IAT) 
         IF (IAT == 0) RETURN  
         CALL HIBFF (IA, IB, IA, IA, 2) 
         IF (IABS(IK1(6)-ID1(6)) > 2) RETURN  
         IF (IABS(IK2(6)-ID2(6)) > 2) RETURN  
         IF (IABS(IK1(5)-ID1(5)) > 4*IK1(3)) RETURN  
         IF (IABS(IK2(5)-ID2(5)) > 4*IK2(3)) RETURN  
         IRS = 0 
         IW1(1,:) = 0 
         IW1(2,:) = 0 
         IW2(1,:) = 0 
         IW2(2,:) = 0 
         IRL = 0 
         LIA2 = LJ(IA)*2 
         LIB2 = LJ(IB)*2 
         IG = MIN0(LIA2,LIB2) + 1 
         DO I = 1, IG, 2 
            KL = I - 1 
            IF (IOCASE == 1) THEN 
! Spin - spin      operator
               CALL SS1212 (IG, KL, IA, IB) 
            ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
               CALL SOO1212 (IG, KL, IA, IB) 
!            CALL SOO1212P(IG,KL,IA,IB)
            ENDIF 
         END DO 
         IP = IABS(LJ(IB)-LJ(IA)) + 1 
         IG = LJ(IB) + LJ(IA) + 1 
         DO I = IP, IG, 2 
            KL = I - 1 
            IF (IOCASE == 1) THEN 
! Spin-spin operator
               CALL SS1221 (IG, KL, IA, IB) 
            ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
               CALL SOOC (IG, KL, IA, IB, IB, IB, IA, IB, IB, IA, KL, TWO13) 
            ENDIF 
         END DO 
      ENDIF 
      RETURN  
      END SUBROUTINE TWO1 
