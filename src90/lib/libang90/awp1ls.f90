!
!     -------------------------------------------------------------
!      A W P 1 L S
!     -------------------------------------------------------------
!
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!                                                                  *
!                N      (ls)  (k1 k3) (k2 k4)  N'     +-           *
!     ELEMENT: (l QLS ::[A  * W      ]      ::l  QLS) -+           *
!                                                     ++           *
!                                                     -- B17 (2.3) *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE AWP1LS(IK, BK, ID, BD, K1, K2, K3, BK4, QM1, QM2, QM3, AW) 
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
      USE awp1g_I 
      USE izas1_I 
      USE itls2_I 
      USE rumt_I 
      USE ixjtik_I 
      USE c0t5s_I 
      USE sls_I 
      USE w1_I 
      USE sixj_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: K1 
      INTEGER  :: K2 
      INTEGER  :: K3 
      REAL(DOUBLE)  :: BK4 
      REAL(DOUBLE)  :: QM1 
      REAL(DOUBLE)  :: QM2 
      REAL(DOUBLE)  :: QM3 
      REAL(DOUBLE)  :: AW 
      INTEGER  :: IK(7) 
      INTEGER  :: ID(7) 
      REAL(DOUBLE)  :: BK(3) 
      REAL(DOUBLE)  :: BD(3) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(7) :: IBT 
      INTEGER :: IQ2, IQ3, IQ, KK1, KK2, KK3, KK4, ITP, ITG, IQM, IT, IE 
      REAL(DOUBLE), DIMENSION(3) :: BT 
      REAL(DOUBLE) :: ENQP, D1, S, W, SI1, SI2 
!-----------------------------------------------
      AW = ZERO 
      IF (ID(3) == 3) THEN 
         IF (MAX0(IK(4),ID(4)) < 3) THEN 
            IF (IK(1) < 300) CALL MES (54) 
            IF (ID(1) < 300) CALL MES (54) 
            CALL AWP1G (K1, K2, K3, BK4, QM1, QM2, QM3, IK, BK, ID, BD, AW) 
            RETURN  
         ENDIF 
      ELSE IF (ID(3) > 3) THEN 
         CALL AWP1G (K1, K2, K3, BK4, QM1, QM2, QM3, IK, BK, ID, BD, AW) 
         RETURN  
      ENDIF 
      IF (IZAS1(ID(7),BD(3),IK(7),BK(3)) == 0) RETURN  
      ENQP = ZERO 
      IQ2 = QM2*TWO + QM2*TENTH 
      IQ3 = QM3*TWO + QM3*TENTH 
      IQ = IQ2 + IQ3 
      KK1 = K1*2 
      KK2 = K2*2 
      KK3 = K3*2 
      KK4 = BK4 + BK4 + TENTH*BK4 
      IF (ITLS2(IK,ID,KK2,KK4,BD,IBT,BT,ITP,ITG,IQ) == 0) RETURN  
      IQM = TWO*DABS(BT(3)) + TENTH 
      DO IT = ITP, ITG 
         CALL RUMT (IT, IBT(3), IBT(7), IBT(6), IBT(5)) 
         IF (IQM > IBT(7)) CYCLE  
         IF (IXJTIK(2*IK(3),KK1,KK2,ID(5),IK(5),IBT(5)) == 0) CYCLE  
         IF (IXJTIK(1,KK3,KK4,ID(6),IK(6),IBT(6)) == 0) CYCLE  
         IBT(1) = IT 
         BT(2) = DBLE(IBT(6))/TWO 
         BT(1) = DBLE(IBT(7))/TWO 
         CALL C0T5S (BT(1), BT(3), QM1, BK(1), BK(3), D1) 
         IF (DABS(D1) < EPS) CYCLE  
         CALL SLS (IK(3), IK(1), IK(7), IK(5), IK(6), IBT(1), IBT(7), IBT(5), &
            IBT(6), S) 
         IF (DABS(S) < EPS) CYCLE  
         CALL W1 (IBT, BT, ID, BD, K1, K3, QM2, QM3, W) 
         IF (DABS(W) < EPS) CYCLE  
         D1 = D1*W*S 
         CALL SIXJ (2*IK(3), KK1, KK2, ID(5), IK(5), IBT(5), 0, SI1) 
         CALL SIXJ (1, KK3, KK4, ID(6), IK(6), IBT(6), 0, SI2) 
         D1 = D1*SI1*SI2/DSQRT(DBLE(IK(7)+1)) 
         ENQP = ENQP + D1 
      END DO 
      AW = ENQP*DSQRT(DBLE((KK2 + 1)*(KK4 + 1))) 
      IE = KK2 + IK(6) + ID(6) + 2 + KK4 + IK(5) + ID(5) 
      IF ((IE/4)*4 /= IE) AW = -AW 
      RETURN  
      END SUBROUTINE AWP1LS 
