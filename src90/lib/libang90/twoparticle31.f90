!
!     -------------------------------------------------------------
!      T W O P A R T I C L E 3 1
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
!                                              N'2 = N2 + 1        *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE TWOPARTICLE31(IA, IB, IIA, IIB, IIC, IID) 
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
      USE two31_I 
      USE hibff_I 
      USE ittk_I 
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
      INTEGER  :: IIA 
      INTEGER  :: IIB 
      INTEGER  :: IIC 
      INTEGER  :: IID 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IP, IKK, IG, I1, I2, I3, I4, I, KL 
!-----------------------------------------------
      IF (IHSH <= 1) RETURN  
      CALL HIBFF (IA, IB, IA, IA, 2) 
      IF (ID1(4) < 2) RETURN  
      IF (ITTK(1,IK2(6),ID2(6)) == 0) RETURN  
      IF (ITTK(2*IK2(3),IK2(5),ID2(5)) == 0) RETURN  
      IF (IABS(IK1(6)-ID1(6)) > 3) RETURN  
      IF (IABS(IK1(5)-ID1(5)) > 6*IK1(3)) RETURN  
      IF (IOCASE == 1) THEN 
! Spin-Spin operator
         IP = ITREXG2(LJ(IA),LJ(IA),LJ(IA),LJ(IB),IKK) + 1 
      ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
         IP = ITREXG(LJ(IA),LJ(IA),LJ(IA),LJ(IB),IKK) + 1 
      ENDIF 
      IF (IKK <= 0) RETURN  
      IG = IP + IKK - 1 
      IRS(:,1) = 0 
      IRL(:,1) = 0 
      IWAA = 0 
      DO I = IP, IG, 2 
         KL = I - 1 
!   COULOMB
!        CALL TWO31(KL,0,KL,0,0,IA,IB,IIA,IIB,IIC,IID,C1,C2)
!   COULOMB
!
         IF (IOCASE == 1) THEN 
! Spin-spin operator
            CALL SSC (IG, KL, IA, IB, IB, IB, IB, IA, IA, IA, KL, TWO31) 
         ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
            CALL SOOC (IG, KL, IA, IB, IB, IB, IB, IA, IA, IA, KL, TWO31) 
         ENDIF 
!
      END DO 
      RETURN  
      END SUBROUTINE TWOPARTICLE31 
