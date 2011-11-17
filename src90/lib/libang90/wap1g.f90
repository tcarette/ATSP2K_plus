!
!     -------------------------------------------------------------
!      W A P 1 G
!     -------------------------------------------------------------
!
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!                                                                  *
!               N      (k1 K2) N'                    +-            *
!     ELEMENT (l QLS::W     ::l QLS)                 -+            *
!                                                    ++            *
!                                                    -- B17 (2.4)  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vilnius,  Lithuania                          December 1993   *
!
      SUBROUTINE WAP1G(K1, K2, K3, BK4, QM1, QM2, QM3, IK, BK, ID, BD, WW) 
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
      USE izas1_I 
      USE itls3_I 
      USE rumt_I 
      USE ixjtik_I 
      USE a1_I 
      USE w1g_I 
      USE sixj_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: K1 
      INTEGER , INTENT(IN) :: K2 
      INTEGER  :: K3 
      REAL(DOUBLE) , INTENT(IN) :: BK4 
      REAL(DOUBLE)  :: QM1 
      REAL(DOUBLE)  :: QM2 
      REAL(DOUBLE)  :: QM3 
      REAL(DOUBLE) , INTENT(OUT) :: WW 
      INTEGER  :: IK(7) 
      INTEGER  :: ID(7) 
      REAL(DOUBLE)  :: BK(3) 
      REAL(DOUBLE)  :: BD(3) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(7) :: IBT 
      INTEGER :: KK1, KK2, KK3, KK4, IQ3, ITP, ITG, IQM, IT, IE 
      REAL(DOUBLE), DIMENSION(3) :: BT 
      REAL(DOUBLE) :: ENQP, D1, W, SI1, SI2 
!-----------------------------------------------
      WW = ZERO 
      IF (IZAS1(ID(7),BD(3),IK(7),BK(3)) == 0) RETURN  
      ENQP = ZERO 
      KK1 = K1*2 
      KK2 = K2*2 
      KK3 = K3*2 
      KK4 = BK4 + BK4 + TENTH*BK4 
      IQ3 = QM3*TWO + QM3*TENTH 
      IF (ID(3) > 9) RETURN  
      IF (ITLS3(IK,ID,KK2,KK4,BK,BD,IBT,BT,ITP,ITG,IQ3) == 0) RETURN  
      IQM = TWO*DABS(BT(3)) + TENTH 
      DO IT = ITP, ITG 
         CALL RUMT (IT, IBT(3), IBT(7), IBT(6), IBT(5)) 
         IF (IQM > IBT(7)) CYCLE  
         IF (IXJTIK(KK1,IK(3)*2,KK2,ID(5),IK(5),IBT(5)) == 0) CYCLE  
         IF (IXJTIK(KK3,1,KK4,ID(6),IK(6),IBT(6)) == 0) CYCLE  
         IBT(1) = IT 
         BT(2) = DBLE(IBT(6))/TWO 
         BT(1) = DBLE(IBT(7))/TWO 
         CALL A1 (IBT, BT, ID, BD, QM3, D1) 
         IF (DABS(D1) < EPS) CYCLE  
         CALL W1G (K1, K3, QM1, QM2, IK, BK, IBT, BT, W) 
         IF (DABS(W) < EPS) CYCLE  
         D1 = D1*W 
         CALL SIXJ (KK1, IK(3)*2, KK2, ID(5), IK(5), IBT(5), 0, SI1) 
         CALL SIXJ (KK3, 1, KK4, ID(6), IK(6), IBT(6), 0, SI2) 
         ENQP = ENQP + D1*SI1*SI2 
      END DO 
      WW = ENQP*DSQRT(DBLE((KK2 + 1)*(KK4 + 1))) 
      IE = KK2 + IK(6) + ID(6) + KK4 + IK(5) + ID(5) 
      IF ((IE/4)*4 /= IE) WW = -WW 
      RETURN  
      END SUBROUTINE WAP1G 
