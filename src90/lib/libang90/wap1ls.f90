!
!     -------------------------------------------------------------
!      W A P 1 L S
!     -------------------------------------------------------------
!
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!                                                                  *
!                  N       (k1)   (ls) (k2)   N'     +-            *
!     ELEMENT:   (l QLS::[ W    * A   ]    ::l QLS)  -+            *
!                                                    ++            *
!                                                    -- B17 (2.3)  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE WAP1LS(IK, BK, ID, BD, K1, K2, K3, BK4, QM1, QM2, QM3, WA) 
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
      USE wap1g_I 
      USE izas1_I 
      USE itls3_I 
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
      REAL(DOUBLE)  :: WA 
      INTEGER  :: IK(7) 
      INTEGER  :: ID(7) 
      REAL(DOUBLE)  :: BK(3) 
      REAL(DOUBLE)  :: BD(3) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(7) :: IBT 
      INTEGER :: KK1, KK2, KK3, IQ3, KK4, ITP, ITG, IQM, IT, IE 
      REAL(DOUBLE), DIMENSION(3) :: BT 
      REAL(DOUBLE) :: ENQP, D1, S, W, SI1, SI2 
!-----------------------------------------------
      WA = ZERO 
      IF (ID(3) == 3) THEN 
         IF (MAX0(IK(4),ID(4)) < 3) THEN 
            IF (IK(1) < 300) CALL MES (55) 
            IF (ID(1) < 300) CALL MES (55) 
            CALL WAP1G (K1, K2, K3, BK4, QM1, QM2, QM3, IK, BK, ID, BD, WA) 
            RETURN  
         ENDIF 
      ELSE IF (ID(3) > 3) THEN 
         CALL WAP1G (K1, K2, K3, BK4, QM1, QM2, QM3, IK, BK, ID, BD, WA) 
         RETURN  
      ENDIF 
      IF (IZAS1(ID(7),BD(3),IK(7),BK(3)) == 0) RETURN  
      ENQP = ZERO 
      KK1 = K1*2 
      KK2 = K2*2 
      KK3 = K3*2 
      IQ3 = QM3*TWO + QM3*TENTH 
      KK4 = BK4 + BK4 + TENTH*BK4 
      IF (ITLS3(IK,ID,KK2,KK4,BK,BD,IBT,BT,ITP,ITG,IQ3) == 0) RETURN  
      IQM = TWO*ABS(BT(3)) + TENTH 
      DO IT = ITP, ITG 
         CALL RUMT (IT, IBT(3), IBT(7), IBT(6), IBT(5)) 
         IF (IQM > IBT(7)) CYCLE  
         IF (IXJTIK(KK1,2*IK(3),KK2,ID(5),IK(5),IBT(5)) == 0) CYCLE  
         IF (IXJTIK(KK3,1,KK4,ID(6),IK(6),IBT(6)) == 0) CYCLE  
         IBT(1) = IT 
         BT(2) = DBLE(IBT(6))/TWO 
         BT(1) = DBLE(IBT(7))/TWO 
         CALL C0T5S (BD(1), BD(3), QM3, BT(1), BT(3), D1) 
         IF (ABS(D1) < EPS) CYCLE  
         CALL SLS (ID(3), IBT(1), IBT(7), IBT(5), IBT(6), ID(1), ID(7), ID(5), &
            ID(6), S) 
         IF (ABS(S) < EPS) CYCLE  
         D1 = D1*S 
         CALL W1 (IK, BK, IBT, BT, K1, K3, QM1, QM2, W) 
         IF (ABS(W) < EPS) CYCLE  
         D1 = D1*W 
         CALL SIXJ (KK1, 2*IK(3), KK2, ID(5), IK(5), IBT(5), 0, SI1) 
         CALL SIXJ (KK3, 1, KK4, ID(6), IK(6), IBT(6), 0, SI2) 
         D1 = D1*SI1*SI2/SQRT(DBLE(IBT(7)+1)) 
         ENQP = ENQP + D1 
      END DO 
      WA = ENQP*SQRT(DBLE((KK2 + 1)*(KK4 + 1))) 
      IE = KK2 + IK(6) + ID(6) + 2 + KK4 + IK(5) + ID(5) 
      IF ((IE/4)*4 /= IE) WA = -WA 
      RETURN  
      END SUBROUTINE WAP1LS 
