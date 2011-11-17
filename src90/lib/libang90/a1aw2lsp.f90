!
!     -------------------------------------------------------------
!      A 1 A W 2 L S P
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES                           *
!                                        N2      (l2,s)  N2'       *
!     OF FOLLOWING MATRIX ELEMENTS:   (nl Q L S::A(2)::nl Q'L'S')* *
!                                        2 2 2 2         2 2 2 2   *
!                                                                  *
!         N1       (l1,s) (k1,k2) (k3,k4) N1'                   +- *
!     *(nl Q L S::[A(1) * W(11)  ]    ::nl  Q'L'S')             -+ *
!         1 1 1 1                         1  1 1 1              ++ *
!                                                               -- *
!                                                                  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE A1AW2LSP(IAA, IBB, KL1, KS1, KL2, BKKS2, QM1, QM2, QM3, QM4, &
         WW) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE MEDEFN_C 
      USE TRK_C 
      USE KAMPAS_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:33:24  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE mes_I 
      USE c0t5s_I 
      USE sls_I 
      USE awp1ls_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IAA 
      INTEGER , INTENT(IN) :: IBB 
      INTEGER  :: KL1 
      INTEGER  :: KS1 
      INTEGER  :: KL2 
      REAL(DOUBLE)  :: BKKS2 
      REAL(DOUBLE)  :: QM1 
      REAL(DOUBLE)  :: QM2 
      REAL(DOUBLE)  :: QM3 
      REAL(DOUBLE)  :: QM4 
      REAL(DOUBLE) , INTENT(OUT) :: WW 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KKS2, IQMM1, IQMM2, IQMM3, IQMM4, IQMM34, NN, IB1, II 
      REAL(DOUBLE) :: A1, S, AW 
!-----------------------------------------------
      WW = ZERO 
      KKS2 = BKKS2 + HALF + TENTH 
      IF (IWAA(KS1+1,KKS2,KL1+1,KL2+1) == 0) THEN 
         RWAA(KS1+1,KKS2,KL1+1,KL2+1) = ZERO 
         IWAA(KS1+1,KKS2,KL1+1,KL2+1) = 1 
         IF (IK1(3) > 3) THEN 
            IF (IK1(4) > 2) CALL MES (33) 
            IF (ID1(4) > 2) CALL MES (33) 
         ENDIF 
         IF (IK2(3) > 3) THEN 
            IF (IK2(4) > 2) CALL MES (33) 
            IF (ID2(4) > 2) CALL MES (33) 
         ENDIF 
         IQMM1 = QM1 + QM1 + TENTH*QM1 
         IF (IK2(4) /= ID2(4) + IQMM1) RETURN  
         IQMM2 = QM2 + QM2 + TENTH*QM2 
         IQMM3 = QM3 + QM3 + TENTH*QM3 
         IQMM4 = QM4 + QM4 + TENTH*QM4 
         IQMM34 = IQMM2 + IQMM3 + IQMM4 
         IF (IK1(4) /= ID1(4) + IQMM34) RETURN  
         CALL C0T5S (BD2(1), BD2(3), QM1, BK2(1), BK2(3), A1) 
         IF (DABS(A1) < EPS) RETURN  
         CALL SLS (IK2(3), IK2(1), IK2(7), IK2(5), IK2(6), ID2(1), ID2(7), ID2(&
            5), ID2(6), S) 
         IF (DABS(S) < EPS) RETURN  
         CALL AWP1LS (IK1, BK1, ID1, BD1, KL1, KL2, KS1, BKKS2, QM2, QM3, QM4, &
            AW) 
         IF (DABS(AW) < EPS) RETURN  
         WW = -A1*AW*S/DSQRT(DBLE(IK2(7)+1)) 
         NN = 1 
         IB1 = IBB - 1 
         NN = SUM(NOSH1(IAA:IB1)) + NN 
         IF (MOD(NN,2) /= 0) WW = -WW 
         RWAA(KS1+1,KKS2,KL1+1,KL2+1) = WW 
      ELSE 
         WW = RWAA(KS1+1,KKS2,KL1+1,KL2+1) 
      ENDIF 
      RETURN  
      END SUBROUTINE A1AW2LSP 
