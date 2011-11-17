!
!     -------------------------------------------------------------
!      W 1 G
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
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE W1G(K1, K2, QM1, QM2, IK, BK, ID, BD, WW) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE RIBOLSF_C 
      USE RIBOLS3_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:11:51  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE izas1_I 
      USE ittk_I 
      USE mes_I 
      USE rumt_I 
      USE ixjtik_I 
      USE a1_I 
      USE sixj_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: K1 
      INTEGER , INTENT(IN) :: K2 
      REAL(DOUBLE)  :: QM1 
      REAL(DOUBLE)  :: QM2 
      REAL(DOUBLE) , INTENT(OUT) :: WW 
      INTEGER  :: IK(7) 
      INTEGER  :: ID(7) 
      REAL(DOUBLE)  :: BK(3) 
      REAL(DOUBLE)  :: BD(3) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(7) :: IBT 
      INTEGER :: KK1, KK2, IQ, ITK, ITD, ITP1, ITP, ITG1, ITG, IQM, IT, IE 
      REAL(DOUBLE), DIMENSION(3) :: BT 
      REAL(DOUBLE) :: ENQP, D1, W, SI1, SI2 
!-----------------------------------------------
      WW = ZERO 
      IF (IZAS1(ID(7),BD(3),IK(7),BK(3)) == 0) RETURN  
      ENQP = ZERO 
      KK1 = K1*2 
      KK2 = K2*2 
      IQ = QM2*TWO + QM2*TENTH 
      IF (ID(3) > 9) RETURN  
      IF (ITTK(ID(5),IK(5),KK1) == 0) RETURN  
      IF (ITTK(ID(6),IK(6),KK2) == 0) RETURN  
      ITK = IK(1) 
      ITD = ID(1) 
      IF (ID(3) == 3) THEN 
         IF (ID(4) > 2) CALL MES (1) 
         IF (IK(4) > 2) CALL MES (1) 
         ITK = ITK - 300 
         ITD = ITD - 300 
         ITP1 = IMPNLSF(ITK) 
         ITP = IMPNLSF(ITD) 
         IF (ITP1 /= ITP) RETURN  
         ITG1 = IMGNLSF(ITK) 
         ITG = IMGNLSF(ITD) 
      ELSE 
         IF (ID(4) > 2) CALL MES (1) 
         IF (IK(4) > 2) CALL MES (1) 
         ITP1 = IMPNLS3(ITK) 
         ITP = IMPNLS3(ITD) 
         IF (ITP1 /= ITP) RETURN  
         ITG1 = IMGNLS3(ITK) 
         ITG = IMGNLS3(ITD) 
      ENDIF 
      IF (ITG1 /= ITG) RETURN  
      IBT(2) = ID(2) 
      IBT(3) = ID(3) 
      IBT(4) = ID(4) + IQ 
      BT(3) = BD(3) + HALF*DBLE(IQ) 
      IQM = TWO*DABS(BT(3)) + TENTH 
      DO IT = ITP, ITG 
         CALL RUMT (IT, IBT(3), IBT(7), IBT(6), IBT(5)) 
         IF (IQM > IBT(7)) CYCLE  
         IF (IXJTIK(IK(3)*2,IK(3)*2,KK1,ID(5),IK(5),IBT(5)) == 0) CYCLE  
         IF (IXJTIK(1,1,KK2,ID(6),IK(6),IBT(6)) == 0) CYCLE  
         IBT(1) = IT 
         BT(2) = DBLE(IBT(6))/TWO 
         BT(1) = DBLE(IBT(7))/TWO 
         CALL A1 (IK, BK, IBT, BT, QM1, D1) 
         IF (DABS(D1) <= EPS) CYCLE  
         CALL A1 (IBT, BT, ID, BD, QM2, W) 
         IF (DABS(W) <= EPS) CYCLE  
         D1 = D1*W 
         CALL SIXJ (IK(3)*2, IK(3)*2, KK1, ID(5), IK(5), IBT(5), 0, SI1) 
         CALL SIXJ (1, 1, KK2, ID(6), IK(6), IBT(6), 0, SI2) 
         ENQP = ENQP + D1*SI1*SI2 
      END DO 
      WW = ENQP*DSQRT(DBLE((KK1 + 1)*(KK2 + 1))) 
      IE = KK1 + IK(6) + ID(6) + KK2 + IK(5) + ID(5) 
      IF ((IE/4)*4 /= IE) WW = -WW 
      RETURN  
      END SUBROUTINE W1G 
