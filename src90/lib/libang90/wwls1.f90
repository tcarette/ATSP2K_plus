!
!     -------------------------------------------------------------
!      W W L S 1
!     -------------------------------------------------------------
!
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!                                                                  *
!                N      (k2 K3) (k2 K3) (0 0) N'     +-            *
!     ELEMENT  (l QLS::[W   *  W    ]      ::l QLS)  -+            *
!                                                    ++            *
!                                                    -- B17 (2.4)  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vilnius,  Lithuania                             March 1995   *
!
      SUBROUTINE WWLS1(IK, BK, ID, BD, K2, K3, QM1, QM2, QM3, QM4, WW) 
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
      USE itls_I 
      USE rumt_I 
      USE ixjtik_I 
      USE w1_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: K2 
      INTEGER  :: K3 
      REAL(DOUBLE)  :: QM1 
      REAL(DOUBLE)  :: QM2 
      REAL(DOUBLE)  :: QM3 
      REAL(DOUBLE)  :: QM4 
      REAL(DOUBLE) , INTENT(OUT) :: WW 
      INTEGER  :: IK(7) 
      INTEGER  :: ID(7) 
      REAL(DOUBLE)  :: BK(3) 
      REAL(DOUBLE)  :: BD(3) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(7) :: IBT 
      INTEGER :: KK2, KK3, IQ3, IQ4, IQ, KK6, KK7, ITP, ITG, IE1, IQM, IT 
      REAL(DOUBLE), DIMENSION(3) :: BT 
      REAL(DOUBLE) :: ENQP, D1, W 
!-----------------------------------------------
      WW = ZERO 
      IF (ID(5) /= IK(5)) RETURN  
      IF (ID(6) /= IK(6)) RETURN  
      IF (IZAS1(ID(7),BD(3),IK(7),BK(3)) == 0) RETURN  
      ENQP = ZERO 
      KK2 = K2*2 
      KK3 = K3*2 
      IQ3 = QM3*TWO + QM3*TENTH 
      IQ4 = QM4*TWO + QM4*TENTH 
      IQ = IQ3 + IQ4 
      IF (ITLS(IK,ID,0,0,BD,IBT,BT,KK6,KK7,ITP,ITG,IQ) == 0) RETURN  
      IE1 = KK2 - IK(5) + KK3 - IK(6) 
      IQM = TWO*DABS(BT(3)) + TENTH 
      DO IT = ITP, ITG 
         CALL RUMT (IT, IBT(3), IBT(7), IBT(6), IBT(5)) 
         IF (IQM > IBT(7)) CYCLE  
         IF (IXJTIK(KK2,KK2,0,ID(5),IK(5),IBT(5)) == 0) CYCLE  
         IF (IXJTIK(KK3,KK3,0,ID(6),IK(6),IBT(6)) == 0) CYCLE  
         IBT(1) = IT 
         BT(2) = DBLE(IBT(6))/TWO 
         BT(1) = DBLE(IBT(7))/TWO 
         CALL W1 (IK, BK, IBT, BT, K2, K3, QM1, QM2, D1) 
         IF (DABS(D1) <= EPS) CYCLE  
         CALL W1 (IBT, BT, ID, BD, K2, K3, QM3, QM4, W) 
         IF (DABS(W) <= EPS) CYCLE  
         D1 = D1*W 
         IF (MOD(IE1 + IBT(5)+IBT(6),4) /= 0) D1 = -D1 
         ENQP = ENQP + D1 
      END DO 
      WW = ENQP/DSQRT(DBLE((KK2 + 1)*(IK(5)+1)*(KK3+1)*(IK(6)+1))) 
      RETURN  
      END SUBROUTINE WWLS1 
