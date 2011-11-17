 
 
!
!     -------------------------------------------------------------
!      W 1 W 2 L S P
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!     ELEMENTS:                                                    *
!                                                                  *
!        N1       (k1)    N1'                                      *
!     (nl Q L S::W(11)::nl Q'L'S') *                               *
!        1 1 1 1          1 1 1 1                                  *
!                                     N2       (k2)    N2'      +- *
!                                * (nl Q L S::W(22)::nl Q'L'S') -+ *
!                                     2 2 2 2          2 2 2 2  ++ *
!                                                               -- *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE W1W2LSP(KL1, KS1, KL2, KS2, QM1, QM2, QM3, QM4, AA) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE TRK_C 
      USE KAMPAS_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:11:51  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE mes_I 
      USE w1_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: KL1 
      INTEGER  :: KS1 
      INTEGER  :: KL2 
      INTEGER  :: KS2 
      REAL(DOUBLE)  :: QM1 
      REAL(DOUBLE)  :: QM2 
      REAL(DOUBLE)  :: QM3 
      REAL(DOUBLE)  :: QM4 
      REAL(DOUBLE) , INTENT(OUT) :: AA 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IQMM1, IQMM2, IQMM3, IQMM4 
      REAL(DOUBLE) :: A1, W 
!-----------------------------------------------
      AA = ZERO 
      IF (IW1(KS1+1,KL1+1) == 0) THEN 
         IW1(KS1+1,KL1+1) = 1 
         IF (IK1(3) > 3) THEN 
            IF (IK1(4) > 2) CALL MES (35) 
            IF (ID1(4) > 2) CALL MES (35) 
         ENDIF 
         IQMM1 = QM1 + QM1 + TENTH*QM1 
         IQMM2 = QM2 + QM2 + TENTH*QM2 
         IF (IK1(4) /= ID1(4) + IQMM1 + IQMM2) THEN 
            A1 = ZERO 
         ELSE 
            CALL W1 (IK1, BK1, ID1, BD1, KL1, KS1, QM1, QM2, A1) 
         ENDIF 
         RW1(KS1+1,KL1+1) = A1 
      ELSE 
         A1 = RW1(KS1+1,KL1+1) 
      ENDIF 
      IF (DABS(A1) < EPS) RETURN  
!
      IF (IW2(KS2+1,KL2+1) == 0) THEN 
         IW2(KS2+1,KL2+1) = 1 
         IF (IK2(3) > 3) THEN 
            IF (IK2(4) > 2) CALL MES (35) 
            IF (ID2(4) > 2) CALL MES (35) 
         ENDIF 
         IQMM3 = QM3 + QM3 + TENTH*QM3 
         IQMM4 = QM4 + QM4 + TENTH*QM4 
         IF (IK2(4) /= ID2(4) + IQMM3 + IQMM4) THEN 
            W = ZERO 
         ELSE 
            CALL W1 (IK2, BK2, ID2, BD2, KL2, KS2, QM3, QM4, W) 
         ENDIF 
         RW2(KS2+1,KL2+1) = W 
      ELSE 
         W = RW2(KS2+1,KL2+1) 
      ENDIF 
      IF (DABS(W) < EPS) RETURN  
      AA = A1*W 
      RETURN  
      END SUBROUTINE W1W2LSP 
