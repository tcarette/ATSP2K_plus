!
!     -------------------------------------------------------------
!      S S C
!     -------------------------------------------------------------
!                                                                  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
!                                                                  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University Nashville, USA           October 1995  *
!
      SUBROUTINE SSC(IG, KL, IA, IB, IC, ID, IIA, IIB, IIC, IID, IREZ, XXX) 
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
      USE savenon_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: IG 
      INTEGER  :: KL 
      INTEGER  :: IA 
      INTEGER  :: IB 
      INTEGER  :: IC 
      INTEGER  :: ID 
      INTEGER , INTENT(IN) :: IIA 
      INTEGER , INTENT(IN) :: IIB 
      INTEGER , INTENT(IN) :: IIC 
      INTEGER , INTENT(IN) :: IID 
      INTEGER  :: IREZ 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KLP, L1, L2, L3, L4, LA, LB, LC, LD 
      REAL(DOUBLE) :: A1, A2, C, CC 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
      KLP = KL + 2 
!      IF(IG-1.LT.KLP) RETURN
      L1 = LJ(IIA) 
      L2 = LJ(IIB) 
      L3 = LJ(IIC) 
      L4 = LJ(IID) 
      CALL SSA (L1, L2, L3, L4, KLP, KL, A1) 
      CALL SSA (L2, L1, L4, L3, KLP, KL, A2) 
      IF (DABS(A1) + DABS(A2) < TWO*EPS) RETURN  
      LA = IJFUL(IIA) 
      LB = IJFUL(IIB) 
      LC = IJFUL(IIC) 
      LD = IJFUL(IID) 
      CALL XXX (KLP, 1, KL, 1, 2, IA, IB, IC, ID, IREZ, C, CC) 
      C = C*A1 
      CC = CC*A2 
      IF (DABS(C) > EPS) CALL SAVENON (8, C, KL, LA, LB, LC, LD, JA, JB, 0) 
      IF (DABS(CC) > EPS) CALL SAVENON (8, CC, KL, LB, LA, LD, LC, JA, JB, 0) 
      RETURN  
      END SUBROUTINE SSC 
