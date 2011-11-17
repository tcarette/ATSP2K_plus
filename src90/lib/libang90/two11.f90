!     -------------------------------------------------------------
!      T W O 1 1
!     -------------------------------------------------------------
!                                                                  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
!                                                                  *
!                                                  1111            *
!                                                                  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE TWO11(KL1, KS1, KL2, KS2, K, IA, C1) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE MEDEFN_C 
      USE DIAGNL_C 
      USE CONSTS_C 
      USE TRK_C 
      USE PERMAT_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:23:13  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ixjtik_I 
      USE w1_I 
      USE sixj_I 
      USE wwpls1_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: KL1 
      INTEGER  :: KS1 
      INTEGER  :: KL2 
      INTEGER  :: KS2 
      INTEGER  :: K 
      INTEGER  :: IA 
      REAL(DOUBLE) , INTENT(OUT) :: C1 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: W, SI1, SI2, WW 
!-----------------------------------------------
      C1 = ZERO 
      W = ZERO 
      IF (IXJTIK(2*KL1,2*KL2,2*K,2*ID1(3),2*ID1(3),2*ID1(3)) /= 0) THEN 
         IF (IXJTIK(2*KS1,2*KS2,2*K,1,1,1) /= 0) THEN 
            CALL W1 (IK1, BK1, ID1, BD1, K, K, HALF, (-HALF), W) 
            IF (DABS(W) > EPS) THEN 
               CALL SIXJ (2*KL1, 2*KL2, 2*K, 2*ID1(3), 2*ID1(3), 2*ID1(3), 0, &
                  SI1) 
               CALL SIXJ (1, 1, 2*K, 2*KS1, 2*KS2, 1, 0, SI2) 
               W = W*SI1*SI2 
            ENDIF 
         ENDIF 
      ENDIF 
      CALL WWPLS1 (KL1, KS1, KL2, KS2, K, K, HALF, (-HALF), HALF, (-HALF), WW) 
      WW = WW/DSQRT(DBLE((2*KL1 + 1)*(2*KL2 + 1)*(2*KS1 + 1)*(2*KS2 + 1))) 
      C1 = HALF*RS(1,1)*RL(1,1)*(WW - W) 
      RETURN  
      END SUBROUTINE TWO11 
