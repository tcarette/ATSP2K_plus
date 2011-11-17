!
!     -------------------------------------------------------------
!      T W O P A R T I C L E 4
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :      N'1 = N1 +- 1        *
!                                             N'2 = N2 +- 1        *
!                                             N'3 = N3 -+ 2        *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!                                                                  *
      SUBROUTINE TWOPARTICLE4(IA, IB, IC, ID) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE MEDEFN_C 
      USE TRK_C 
      USE TRK2_C 
      USE PERMAT_C 
      USE KAMPAS_C 
      USE CASEOP_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:23:13  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE two41_I 
      USE two42_I 
      USE eile_I 
      USE hibff_I 
      USE ittk_I 
      USE rlsp00_I 
      USE itrexg2_I 
      USE itrexg_I 
      USE ssc_I 
      USE sooc_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: IA 
      INTEGER  :: IB 
      INTEGER  :: IC 
      INTEGER  :: ID 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KGAL, IAA, IBB, ICC, IAT, IP, IKK, IG, I1, I2, I, KL 
!-----------------------------------------------
      IF (IHSH <= 2) RETURN  
!
      IF (IOCASE == 1) THEN 
! Spin - spin      operator
         KGAL = 2 
      ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
         KGAL = 1 
      ENDIF 
!
      IF (IA == IB) THEN 
         CALL EILE (IA, IC, ID, IAA, IBB, ICC) 
         CALL HIBFF (IC, ID, IA, IC, 3) 
      ELSE IF (IC == ID) THEN 
         CALL EILE (IA, IB, IC, IAA, IBB, ICC) 
         CALL HIBFF (IA, IB, IC, IA, 3) 
         IF (ID3(4) < 2) RETURN  
      ELSE 
         WRITE (6, '(A)') ' ERROR IN SUBROUTINE TWOPARTICLE4  ' 
         STOP  
      ENDIF 
      IF (ITTK(1,IK1(6),ID1(6)) == 0) RETURN  
      IF (ITTK(1,IK2(6),ID2(6)) == 0) RETURN  
      IF (IABS(IK3(6)-ID3(6)) > 2) RETURN  
      IF (ITTK(2*IK1(3),IK1(5),ID1(5)) == 0) RETURN  
      IF (ITTK(2*IK2(3),IK2(5),ID2(5)) == 0) RETURN  
      IF (IABS(IK3(5)-ID3(5)) > 4*IK3(3)) RETURN  
      CALL RLSP00 (1, IAA, IBB, ICC, ICC, 0, IAT) 
      IF (IAT == 0) RETURN  
      CALL RLSP00 (2, IAA, IBB, ICC, ICC, 2*KGAL, IAT) 
      IF (IAT == 0) RETURN  
      CALL RLSP00 (3, IAA, IBB, ICC, ICC, 2*KGAL, IAT) 
      IF (IAT == 0) RETURN  
      IF (IOCASE == 1) THEN 
! Spin-Spin operator
         IP = ITREXG2(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK) + 1 
!        IP=ITREXG2(LJ(IA),LJ(IC),LJ(IB),LJ(IC),IKK)+1
      ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
         IP = ITREXG(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK) + 1 
!        IP=ITREXG(LJ(IA),LJ(IC),LJ(IB),LJ(IC),IKK)+1
      ENDIF 
      IF (IKK <= 0) RETURN  
      IG = IP + IKK - 1 
      IRS = 0 
      IW1(1,:) = 0 
      IW1(2,:) = 0 
      IRL = 0 
      DO I = IP, IG 
         KL = I - 1 
         IF (IA == IB) THEN 
!
            IF (IOCASE == 1) THEN 
! Spin-spin operator
               CALL SSC (IG, KL, IC, ID, IA, IA, IA, IB, IC, ID, KL, TWO41) 
            ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
               CALL SOOC (IG, KL, IC, ID, IA, IA, IA, IB, IC, ID, KL, TWO41) 
            ENDIF 
!
         ELSE IF (IC == ID) THEN 
!
            IF (IOCASE == 1) THEN 
! Spin-spin operator
               CALL SSC (IG, KL, IA, IB, IC, IC, IA, IB, IC, ID, KL, TWO42) 
            ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
               CALL SOOC (IG, KL, IA, IB, IC, IC, IA, IB, IC, ID, KL, TWO42) 
            ENDIF 
!
         ELSE 
            WRITE (6, '(A)') ' ERROR IN SUBROUTINE TWOPARTICLE4  ' 
            STOP  
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE TWOPARTICLE4 
