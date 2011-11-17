!
!     -------------------------------------------------------------
!      A 1 A W 2 L S
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES                           *
!                                        N1      (l1,s)  N1'       *
!     OF FOLLOWING MATRIX ELEMENTS:   (nl Q L S::A(1)::nl Q'L'S')* *
!                                        1 1 1 1         1 1 1 1   *
!                                                                  *
!         N2       (l2,s) (k1,k2) (l1,s) N2'                    +- *
!     *(nl Q L S::[A(2) * W(22)  ]   ::nl  Q'L'S')              -+ *
!         2 2 2 2                        2  2 2 2               ++ *
!                                                               -- *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!
      SUBROUTINE A1AW2LS(IK1, IK2, BK1, BK2, ID1, ID2, BD1, BD2, K1, K2, QM1, &
         QM2, QM3, QM4, WW) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
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
      INTEGER :: IQMM1, IQMM2, IQMM3, IQMM4, IQMM34 
      REAL(DOUBLE) :: A1, S, BKK4, AW 
!-----------------------------------------------
      WW = ZERO 
      IF (IK1(3) > 3) THEN 
         IF (IK1(4) > 2) CALL MES (33) 
         IF (ID1(4) > 2) CALL MES (33) 
      ENDIF 
      IF (IK2(3) > 3) THEN 
         IF (IK2(4) > 2) CALL MES (33) 
         IF (ID2(4) > 2) CALL MES (33) 
      ENDIF 
      IQMM1 = QM1 + QM1 + TENTH*QM1 
      IF (IK1(4) /= ID1(4) + IQMM1) RETURN  
      IQMM2 = QM2 + QM2 + TENTH*QM2 
      IQMM3 = QM3 + QM3 + TENTH*QM3 
      IQMM4 = QM4 + QM4 + TENTH*QM4 
      IQMM34 = IQMM2 + IQMM3 + IQMM4 
      IF (IK2(4) /= ID2(4) + IQMM34) RETURN  
      CALL C0T5S (BD1(1), BD1(3), QM1, BK1(1), BK1(3), A1) 
      IF (ABS(A1) < EPS) RETURN  
      CALL SLS (IK1(3), IK1(1), IK1(7), IK1(5), IK1(6), ID1(1), ID1(7), ID1(5)&
         , ID1(6), S) 
      IF (ABS(S) < EPS) RETURN  
      BKK4 = HALF 
      CALL AWP1LS (IK2, BK2, ID2, BD2, K1, IK1(3), K2, BKK4, QM2, QM3, QM4, AW) 
      IF (ABS(AW) < EPS) RETURN  
      WW = -A1*AW*S/SQRT(DBLE(IK1(7)+1)) 
      RETURN  
      END SUBROUTINE A1AW2LS 
