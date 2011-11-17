!
!     -------------------------------------------------------------
!      S S 1 2 2 1
!     -------------------------------------------------------------
!                                                                  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
!                                                                  *
!                                                                  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
!
!GG elektrostatine
!      CALL COULOMBLS(LJ(IA),LJ(IB),LJ(IB),LJ(IA),KL1,A1)
!      A1=A1*TWO*DSQRT(DBLE(2*KL1+1))
!GG elektrostatine
      SUBROUTINE SS1221(IG, KL, IA, IB) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE MEDEFN_C 
      USE DIAGNL_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:25:04  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ssa_I 
      USE two13_I 
      USE savenon_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IG 
      INTEGER  :: KL 
      INTEGER  :: IA 
      INTEGER  :: IB 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KLP, LA, LB, INE3 
      REAL(DOUBLE) :: A1, C, CC 
!-----------------------------------------------
      KLP = KL + 2 
      IF (IG - 1 < KLP) RETURN  
      LA = IJFUL(IA) 
      LB = IJFUL(IB) 
      CALL SSA (LJ(IA), LJ(IB), LJ(IB), LJ(IA), KLP, KL, A1) 
      IF (DABS(A1) < EPS) RETURN  
      INE3 = 0 
      CALL TWO13 (KLP, 1, KL, 1, 2, IA, IB, IB, IA, INE3, C, CC) 
      C = C*A1*TWO 
      IF (DABS(C) > EPS) CALL SAVENON (8, C, KL, LA, LB, LB, LA, JA, JB, 0) 
      RETURN  
      END SUBROUTINE SS1221 
