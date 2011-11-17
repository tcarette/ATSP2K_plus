!
!     -------------------------------------------------------------
!      W A 1 A 2 L S
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!                                                                  *
!                        N1       (k1,k2) (l1,s) (l2,s) N1'        *
!     ELEMENT         (nl Q L S::[W(11) * A(1)  ]   ::nl  Q'L'S')* *
!                         1 1 1 1                        1  1 1 1  *
!                                                                  *
!         N2      (l2,s)  N2'                                   +- *
!     *(nl Q L S::A(2)::nl Q'L'S')                              -+ *
!         2 2 2 2         2 2 2 2                               ++ *
!                                                               -- *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!
      SUBROUTINE WA1A2LS(IK1, IK2, BK1, BK2, ID1, ID2, BD1, BD2, K1, K2, QM1, &
         QM2, QM3, QM4, WW) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
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
      INTEGER  :: K1 
      INTEGER  :: K2 
      REAL(DOUBLE)  :: QM1 
      REAL(DOUBLE)  :: QM2 
      REAL(DOUBLE)  :: QM3 
      REAL(DOUBLE)  :: QM4 
      REAL(DOUBLE) , INTENT(OUT) :: WW 
      INTEGER  :: IK1(7) 
      INTEGER  :: IK2(7) 
      INTEGER  :: ID1(7) 
      INTEGER  :: ID2(7) 
      REAL(DOUBLE)  :: BK1(3) 
      REAL(DOUBLE)  :: BK2(3) 
      REAL(DOUBLE)  :: BD1(3) 
      REAL(DOUBLE)  :: BD2(3) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IQMM1, IQMM2, IQMM3, IQMM23, IQMM4 
      REAL(DOUBLE) :: A1, S, BK, WA 
!-----------------------------------------------
      WW = ZERO 
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
      IF (IK1(4) /= ID1(4) + IQMM23) RETURN  
      IQMM4 = QM4 + QM4 + TENTH*QM4 
      IF (IK2(4) /= ID2(4) + IQMM4) RETURN  
      CALL C0T5S (BD2(1), BD2(3), QM4, BK2(1), BK2(3), A1) 
      IF (ABS(A1) < EPS) RETURN  
      CALL SLS (IK2(3), IK2(1), IK2(7), IK2(5), IK2(6), ID2(1), ID2(7), ID2(5)&
         , ID2(6), S) 
      IF (ABS(S) < EPS) RETURN  
      BK = HALF 
      CALL WAP1LS (IK1, BK1, ID1, BD1, K1, IK2(3), K2, BK, QM1, QM2, QM3, WA) 
      IF (ABS(WA) < EPS) RETURN  
      WW = -A1*WA*S/SQRT(DBLE(IK2(7)+1)) 
      RETURN  
      END SUBROUTINE WA1A2LS 
