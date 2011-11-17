!
!     -------------------------------------------------------------
!      W W P L S 1
!     -------------------------------------------------------------
!
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!                                                                  *
!               N      (k1 K2) (k3 K4) (k5 k6) N'    +-            *
!     ELEMENT (l QLS::[W   *  W    ]        ::l QLS) -+            *
!                                                    ++            *
!                                                    -- B17 (2.4)  *
!                                                                  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE WWPLS1(K1, K2, K3, K4, K5, K6, QM1, QM2, QM3, QM4, WW) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE TRK_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:11:51  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE izas1_I 
      USE itls_I 
      USE rumt_I 
      USE ixjtik_I 
      USE w1_I 
      USE sixj_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: K1 
      INTEGER  :: K2 
      INTEGER  :: K3 
      INTEGER  :: K4 
      INTEGER  :: K5 
      INTEGER  :: K6 
      REAL(DOUBLE)  :: QM1 
      REAL(DOUBLE)  :: QM2 
      REAL(DOUBLE)  :: QM3 
      REAL(DOUBLE)  :: QM4 
      REAL(DOUBLE) , INTENT(OUT) :: WW 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(7) :: IBT 
      INTEGER :: KK1, KK2, KK3, KK4, IQ3, IQ4, IQ, KK5, KK6, ITP, ITG, IQM, IT 
      REAL(DOUBLE), DIMENSION(3) :: BT 
      REAL(DOUBLE) :: ENQP, D1, W, SI1, SI2 
!-----------------------------------------------
      WW = ZERO 
      IF (IZAS1(ID1(7),BD1(3),IK1(7),BK1(3)) == 0) RETURN  
      ENQP = ZERO 
      KK1 = K1*2 
      KK2 = K2*2 
      KK3 = K3*2 
      KK4 = K4*2 
      IQ3 = QM3*TWO + QM3*TENTH 
      IQ4 = QM4*TWO + QM4*TENTH 
      IQ = IQ3 + IQ4 
      IF (ITLS(IK1,ID1,K5,K6,BD1,IBT,BT,KK5,KK6,ITP,ITG,IQ) == 0) RETURN  
      IQM = TWO*DABS(BT(3)) + TENTH 
      DO IT = ITP, ITG 
         CALL RUMT (IT, IBT(3), IBT(7), IBT(6), IBT(5)) 
         IF (IQM > IBT(7)) CYCLE  
         IF (IXJTIK(KK1,KK3,KK5,ID1(5),IK1(5),IBT(5)) == 0) CYCLE  
         IF (IXJTIK(KK2,KK4,KK6,ID1(6),IK1(6),IBT(6)) == 0) CYCLE  
         IBT(1) = IT 
         BT(2) = DBLE(IBT(6))/TWO 
         BT(1) = DBLE(IBT(7))/TWO 
         CALL W1 (IK1, BK1, IBT, BT, K1, K2, QM1, QM2, D1) 
         IF (DABS(D1) <= EPS) CYCLE  
         CALL W1 (IBT, BT, ID1, BD1, K3, K4, QM3, QM4, W) 
         IF (DABS(W) <= EPS) CYCLE  
         D1 = D1*W 
         CALL SIXJ (KK1, KK3, KK5, ID1(5), IK1(5), IBT(5), 0, SI1) 
         CALL SIXJ (KK2, KK4, KK6, ID1(6), IK1(6), IBT(6), 0, SI2) 
         ENQP = ENQP + D1*SI1*SI2 
      END DO 
      WW = ENQP*DSQRT(DBLE((KK5 + 1)*(KK6 + 1))) 
      IF (MOD(KK5 + IK1(6)+ID1(6)+KK6+IK1(5)+ID1(5),4) /= 0) WW = -WW 
      RETURN  
      END SUBROUTINE WWPLS1 
