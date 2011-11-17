!
!     -------------------------------------------------------------
!      T W O P A R T I C L E 3 3
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
!                                              N'2 = N2 + 1        *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE TWOPARTICLE33(IA, IB, IC, IIA, IIB, IIC, IID) 
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
      USE two33a_I 
      USE two33b_I 
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
      INTEGER  :: IC 
      INTEGER  :: IIA 
      INTEGER  :: IIB 
      INTEGER  :: IIC 
      INTEGER  :: IID 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I1, I2, IP, IKK, IG, I, KL 
!-----------------------------------------------
      IF (IHSH <= 2) RETURN  
      CALL HIBFF (IA, IB, IC, IA, 3) 
      IF (ITTK(1,IK1(6),ID1(6)) == 0) RETURN  
      IF (ITTK(1,IK2(6),ID2(6)) == 0) RETURN  
      IF (IABS(IK3(6)-ID3(6)) > 2) RETURN  
      IF (ITTK(2*IK1(3),IK1(5),ID1(5)) == 0) RETURN  
      IF (ITTK(2*IK2(3),IK2(5),ID2(5)) == 0) RETURN  
      IF (IABS(IK3(5)-ID3(5)) > 4*IK3(3)) RETURN  
      IRS = 0 
      IW1(1,:) = 0 
      IW1(2,:) = 0 
      IRL = 0 
      IF (IOCASE == 1) THEN 
! Spin-Spin operator
         IP = ITREXG2(LJ(IA),LJ(IB),LJ(IC),LJ(IC),IKK) + 1 
      ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
         IP = ITREXG(LJ(IA),LJ(IB),LJ(IC),LJ(IC),IKK) + 1 
      ENDIF 
      IF (IKK > 0) THEN 
         IG = IP + IKK - 1 
         DO I = IP, IG 
            KL = I - 1 
!
            IF (IOCASE == 1) THEN 
! Spin-spin operator
               CALL SSC (IG, KL, IA, IB, IC, IC, IB, IC, IA, IC, KL, TWO33A) 
            ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
               CALL SOOC (IG, KL, IA, IB, IC, IC, IB, IC, IA, IC, KL, TWO33A) 
            ENDIF 
!
         END DO 
      ENDIF 
      IF (IOCASE == 1) THEN 
! Spin-Spin operator
         IP = ITREXG2(LJ(IA),LJ(IC),LJ(IB),LJ(IC),IKK) + 1 
      ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
         IP = ITREXG(LJ(IA),LJ(IC),LJ(IB),LJ(IC),IKK) + 1 
      ENDIF 
      IF (IKK <= 0) RETURN  
      IG = IP + IKK - 1 
      DO I = IP, IG 
         KL = I - 1 
         IF (IOCASE == 1) THEN 
! Spin-spin operator
            CALL SSC (IG, KL, IA, IB, IC, IC, IC, IB, IA, IC, KL, TWO33B) 
         ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
            CALL SOOC (IG, KL, IA, IB, IC, IC, IC, IB, IA, IC, KL, TWO33B) 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE TWOPARTICLE33 
