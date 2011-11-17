!
!     -------------------------------------------------------------
!      W 1
!     -------------------------------------------------------------
!
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!                                                                  *
!                     N      (K2 K3)  N'      +-                   *
!     ELEMENT:      (l QLS::W      ::l QLS)   -+                   *
!                                             ++                   *
!                                             -- S5(1.47),(1.48),  *
!                                                  (1.49),(1.50).  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE W1(IK, BK, ID, BD, K2, K3, QM1, QM2, W) 
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
      USE w1g_I 
      USE c1e1sm_I 
      USE rwls_I 
      USE cle0sm_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: K2 
      INTEGER  :: K3 
      REAL(DOUBLE)  :: QM1 
      REAL(DOUBLE)  :: QM2 
      REAL(DOUBLE)  :: W 
      INTEGER  :: IK(7) 
      INTEGER  :: ID(7) 
      REAL(DOUBLE)  :: BK(3) 
      REAL(DOUBLE)  :: BD(3) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IQM2, K, IQ, K1 
      REAL(DOUBLE) :: QQ, A, AA, WK1, SS, S2 
!-----------------------------------------------
      W = ZERO 
      IF (ID(3) == 3) THEN 
         IF (MAX0(IK(4),ID(4)) < 3) THEN 
            IF (IK(1) < 300) CALL MES (56) 
            IF (ID(1) < 300) CALL MES (56) 
            IQM2 = QM2 + QM2 + QM2*EPS 
            IF (ID(4) + IQM2 > 2) CALL MES (2) 
            CALL W1G (K2, K3, QM1, QM2, IK, BK, ID, BD, W) 
            RETURN  
         ENDIF 
      ELSE IF (ID(3) > 3) THEN 
         IQM2 = QM2 + QM2 + QM2*EPS 
         IF (ID(4) + IQM2 > 2) CALL MES (2) 
         CALL W1G (K2, K3, QM1, QM2, IK, BK, ID, BD, W) 
         RETURN  
      ENDIF 
      K = K2 + K3 
      QQ = QM1 + QM2 
      IF (ABS(QQ) >= EPS) THEN 
         IF (MOD(K,2) /= 0) RETURN  
         IQ = QQ + QQ*TENTH 
         IF (IK(4) - ID(4) - 2*IQ /= 0) RETURN  
         CALL C1E1SM (BD(1), BD(3), QQ, BK(1), BK(3), A) 
         IF (ABS(A) < EPS) RETURN  
         CALL RWLS (1, K2, K3, IK(3), IK(1), ID(1), AA) 
         A = AA*A 
         W = A/SQRT(TWO*BK(1)+ONE) 
      ELSE 
         IF (IK(4) /= ID(4)) RETURN  
         IF (K /= 0) THEN 
            K1 = 1 
            IF (MOD(K,2) /= 0) K1 = 0 
            WK1 = DBLE(K1) 
            CALL CLE0SM (BD(1), BD(3), WK1, BK(1), BK(3), A) 
            IF (ABS(A) < EPS) RETURN  
            CALL RWLS (K1, K2, K3, IK(3), IK(1), ID(1), AA) 
            A = AA*A 
            W = A/SQRT(FOUR*BK(1)+TWO) 
            IF (QM1 >= EPS) RETURN  
            IF (MOD(K,2) /= 0) W = -W 
         ELSE 
            IF (ID(1) /= IK(1)) RETURN  
            IF (QM1 >= EPS) THEN 
               A = -DBLE(ID(4)) 
            ELSE 
               A = DBLE(4*ID(3)+2-ID(4)) 
            ENDIF 
            SS = DBLE((IK(5)+1)*(IK(6)+1)) 
            S2 = DBLE(4*IK(3)+2) 
            W = A*SQRT(SS/S2) 
         ENDIF 
      ENDIF 
      RETURN  
      END SUBROUTINE W1 
