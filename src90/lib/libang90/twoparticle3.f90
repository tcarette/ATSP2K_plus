!
!     -------------------------------------------------------------
!      T W O P A R T I C L E 3
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
!                                              N'2 = N2 + 1        *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE TWOPARTICLE3(IA, IB, IC, ID) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE MEDEFN_C 
      USE CASEOP_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:23:13  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rlsp0_I 
      USE twoparticle31_I 
      USE twoparticle32_I 
      USE eile_I 
      USE rlsp00_I 
      USE twoparticle33_I 
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
      INTEGER :: KGAL, IIA, IIB, IAT, IAA, IBB, ICC 
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
      IF (IB == ID) THEN 
         IF (IA==IB .OR. IC==IB) THEN 
            IF (IA == IC) GO TO 10 
            IIA = MIN0(IA,IC) 
            IIB = MAX0(IA,IC) 
            CALL RLSP0 (1, IIA, IIB, 0, IAT) 
            IF (IAT == 0) RETURN  
            CALL RLSP0 (2, IIA, IIB, 2*KGAL, IAT) 
            IF (IAT == 0) RETURN  
            CALL RLSP0 (3, IIA, IIB, 2*KGAL, IAT) 
            IF (IAT == 0) RETURN  
            IF (IC == IB) THEN 
               CALL TWOPARTICLE31 (IC, IA, IA, IB, IC, ID) 
            ELSE 
               CALL TWOPARTICLE32 (IC, IA, IA, IB, IC, ID) 
            ENDIF 
         ELSE 
            CALL EILE (IA, IB, IC, IAA, IBB, ICC) 
            CALL RLSP00 (1, IAA, IBB, ICC, ICC, 0, IAT) 
            IF (IAT == 0) RETURN  
            CALL RLSP00 (2, IAA, IBB, ICC, ICC, 2*KGAL, IAT) 
            IF (IAT == 0) RETURN  
            CALL RLSP00 (3, IAA, IBB, ICC, ICC, 2*KGAL, IAT) 
            IF (IAT == 0) RETURN  
            CALL TWOPARTICLE33 (IC, IA, IB, IA, IB, IC, ID) 
         ENDIF 
      ELSE IF (IA == IC) THEN 
         IF (IB==IA .OR. ID==IA) THEN 
            IF (IB == ID) GO TO 10 
            IIA = MIN0(IB,ID) 
            IIB = MAX0(IB,ID) 
            CALL RLSP0 (1, IIA, IIB, 0, IAT) 
            IF (IAT == 0) RETURN  
            CALL RLSP0 (2, IIA, IIB, 2*KGAL, IAT) 
            IF (IAT == 0) RETURN  
            CALL RLSP0 (3, IIA, IIB, 2*KGAL, IAT) 
            IF (IAT == 0) RETURN  
            IF (ID == IA) THEN 
               CALL TWOPARTICLE31 (ID, IB, IA, IB, IC, ID) 
            ELSE 
               CALL TWOPARTICLE32 (ID, IB, IA, IB, IC, ID) 
            ENDIF 
         ELSE 
            CALL EILE (IA, IB, ID, IAA, IBB, ICC) 
            CALL RLSP00 (1, IAA, IBB, ICC, ICC, 0, IAT) 
            IF (IAT == 0) RETURN  
            CALL RLSP00 (2, IAA, IBB, ICC, ICC, 2*KGAL, IAT) 
            IF (IAT == 0) RETURN  
            CALL RLSP00 (3, IAA, IBB, ICC, ICC, 2*KGAL, IAT) 
            IF (IAT == 0) RETURN  
            CALL TWOPARTICLE33 (ID, IB, IA, IA, IB, IC, ID) 
         ENDIF 
      ELSE IF (IA == ID) THEN 
         IF (IB==IA .OR. IC==IA) THEN 
            IF (IB == IC) GO TO 10 
            IIA = MIN0(IB,IC) 
            IIB = MAX0(IB,IC) 
            CALL RLSP0 (1, IIA, IIB, 0, IAT) 
            IF (IAT == 0) RETURN  
            CALL RLSP0 (2, IIA, IIB, 2*KGAL, IAT) 
            IF (IAT == 0) RETURN  
            CALL RLSP0 (3, IIA, IIB, 2*KGAL, IAT) 
            IF (IAT == 0) RETURN  
            IF (IC == ID) THEN 
               CALL TWOPARTICLE31 (IC, IB, IA, IB, IC, ID) 
            ELSE 
               CALL TWOPARTICLE32 (IC, IB, IA, IB, IC, ID) 
            ENDIF 
         ELSE 
            CALL EILE (IA, IB, IC, IAA, IBB, ICC) 
            CALL RLSP00 (1, IAA, IBB, ICC, ICC, 0, IAT) 
            IF (IAT == 0) RETURN  
            CALL RLSP00 (2, IAA, IBB, ICC, ICC, 2*KGAL, IAT) 
            IF (IAT == 0) RETURN  
            CALL RLSP00 (3, IAA, IBB, ICC, ICC, 2*KGAL, IAT) 
            IF (IAT == 0) RETURN  
            CALL TWOPARTICLE33 (IC, IB, IA, IA, IB, ID, IC) 
         ENDIF 
      ELSE IF (IB == IC) THEN 
         IF (IA==IB .OR. ID==IB) THEN 
            IF (IA == ID) GO TO 10 
            IIA = MIN0(IA,ID) 
            IIB = MAX0(IA,ID) 
            CALL RLSP0 (1, IIA, IIB, 0, IAT) 
            IF (IAT == 0) RETURN  
            CALL RLSP0 (2, IIA, IIB, 2*KGAL, IAT) 
            IF (IAT == 0) RETURN  
            CALL RLSP0 (3, IIA, IIB, 2*KGAL, IAT) 
            IF (IAT == 0) RETURN  
            IF (ID == IB) THEN 
               CALL TWOPARTICLE31 (ID, IA, IA, IB, IC, ID) 
            ELSE 
               CALL TWOPARTICLE32 (ID, IA, IA, IB, IC, ID) 
            ENDIF 
         ELSE 
            CALL EILE (IA, IB, ID, IAA, IBB, ICC) 
            CALL RLSP00 (1, IAA, IBB, ICC, ICC, 0, IAT) 
            IF (IAT == 0) RETURN  
            CALL RLSP00 (2, IAA, IBB, ICC, ICC, 2*KGAL, IAT) 
            IF (IAT == 0) RETURN  
            CALL RLSP00 (3, IAA, IBB, ICC, ICC, 2*KGAL, IAT) 
            IF (IAT == 0) RETURN  
            CALL TWOPARTICLE33 (ID, IA, IB, IA, IB, ID, IC) 
         ENDIF 
      ENDIF 
      RETURN  
   10 CONTINUE 
      WRITE (6, '(A)') ' ERROR IN NONRELAT3' 
      STOP  
      END SUBROUTINE TWOPARTICLE3 
