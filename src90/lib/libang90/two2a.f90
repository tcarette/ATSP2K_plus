!
!     -------------------------------------------------------------
!      T W O 2 A
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
      SUBROUTINE TWO2A(IA, IB) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
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
      USE two2_I 
      USE rlsp0_I 
      USE hibff_I 
      USE ss1122_I 
      USE soo1122p_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: IA 
      INTEGER  :: IB 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KGAL, IIA, IIB, IAT, I1, I2, IP, IG, I, KL 
!-----------------------------------------------
      IF (IHSH <= 1) RETURN  
!
      IF (IOCASE == 1) THEN 
! Spin - spin      operator
         KGAL = 2 
      ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
         KGAL = 1 
      ENDIF 
!
      IIA = MIN0(IA,IB) 
      IIB = MAX0(IA,IB) 
      CALL RLSP0 (1, IIA, IIB, 0, IAT) 
      IF (IAT == 0) RETURN  
      CALL RLSP0 (2, IIA, IIB, 2*KGAL, IAT) 
      IF (IAT == 0) RETURN  
      CALL RLSP0 (3, IIA, IIB, 2*KGAL, IAT) 
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
      IP = IABS(LJ(IB)-LJ(IA)) + 1 
      IG = LJ(IB) + LJ(IA) + 1 
      DO I = IP, IG, 2 
         KL = I - 1 
         IF (IOCASE == 1) THEN 
! Spin-spin operator
            CALL SS1122 (IG, KL, IA, IB) 
         ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
            CALL SOO1122P (IG, KL, IA, IB, IB, IB, IA, IA, IB, IB, KL, TWO2) 
!          CALL SOO1122(IG,KL,IA,IB,TWO2)
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE TWO2A 
