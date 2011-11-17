!
!     -------------------------------------------------------------
!      W A 1 A 2 L S P
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!                                                                  *
!                       N2       (k1,k2) (l2,s) (k3,k4) N2'        *
!     ELEMENT        (nl Q L S::[W(22) * A(2)  ]    ::nl  Q'L'S')* *
!                       2 2 2 2                         2  2 2 2   *
!                                                                  *
!         N1      (l1,s)  N1'                                   +- *
!     *(nl Q L S::A(1)::nl Q'L'S')                              -+ *
!         1 1 1 1         1 1 1 1                               ++ *
!                                                               -- *
!                                                                  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE WA1A2LSP(IAA, IBB, KL1, KS1, KL2, BKKS2, QM1, QM2, QM3, QM4, &
         WW) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE MEDEFN_C 
      USE TRK_C 
      USE KAMPAS_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:11:51  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE mes_I 
      USE c0t5s_I 
      USE sls_I 
      USE wap1ls_I 
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
      INTEGER :: KKS2, IQMM1, IQMM2, IQMM3, IQMM23, IQMM4, NN, IB1, II 
      REAL(DOUBLE) :: A1, S, WA 
!-----------------------------------------------
      WW = ZERO 
      KKS2 = BKKS2 + HALF + TENTH 
      IF (IWAA(KS1+1,KKS2,KL1+1,KL2+1) == 0) THEN 
         RWAA(KS1+1,KKS2,KL1+1,KL2+1) = ZERO 
         IWAA(KS1+1,KKS2,KL1+1,KL2+1) = 1 
         IF (IK1(3) > 3) THEN 
            IF (IK1(4) > 2) CALL MES (34) 
            IF (ID1(4) > 2) CALL MES (34) 
         ENDIF 
         IF (IK2(3) > 3) THEN 
            IF (IK2(4) > 2) CALL MES (34) 
            IF (ID2(4) > 2) CALL MES (34) 
         ENDIF 
         IQMM1 = QM1 + QM1 + TENTH*QM1 
         IQMM2 = QM2 + QM2 + TENTH*QM2 
         IQMM3 = QM3 + QM3 + TENTH*QM3 
         IQMM23 = IQMM1 + IQMM2 + IQMM3 
         IF (IK2(4) /= ID2(4) + IQMM23) RETURN  
         IQMM4 = QM4 + QM4 + TENTH*QM4 
         IF (IK1(4) /= ID1(4) + IQMM4) RETURN  
         CALL C0T5S (BD1(1), BD1(3), QM4, BK1(1), BK1(3), A1) 
         IF (DABS(A1) < EPS) RETURN  
         CALL SLS (IK1(3), IK1(1), IK1(7), IK1(5), IK1(6), ID1(1), ID1(7), ID1(&
            5), ID1(6), S) 
         IF (DABS(S) < EPS) RETURN  
         CALL WAP1LS (IK2, BK2, ID2, BD2, KL1, KL2, KS1, BKKS2, QM1, QM2, QM3, &
            WA) 
         IF (DABS(WA) < EPS) RETURN  
         WW = -A1*WA*S/DSQRT(DBLE(IK1(7)+1)) 
         NN = 1 
         IB1 = IBB - 1 
         NN = SUM(NOSH1(IAA:IB1)) + NN 
         IF (MOD(NN,2) /= 0) WW = -WW 
         RWAA(KS1+1,KKS2,KL1+1,KL2+1) = WW 
      ELSE 
         WW = RWAA(KS1+1,KKS2,KL1+1,KL2+1) 
      ENDIF 
      RETURN  
      END SUBROUTINE WA1A2LSP 
