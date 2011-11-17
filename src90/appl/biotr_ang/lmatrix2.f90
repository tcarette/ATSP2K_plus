!
!     -------------------------------------------------------------
!      L M A T R I X 2
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF ONE PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
!                                              N'2 = N2 + 1        *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville           Sebtember 1997   *
!
      SUBROUTINE LMATRIX2(IA, IB) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE MEDEFN_C 
      USE DIAGNL_C 
      USE TRK_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:58:24  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE recoupls0_I 
      USE recoupls2_I 
      USE hibff_I 
      USE a1a2ls_I 
      USE savenon_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: IA 
      INTEGER  :: IB 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IAA, IBB, IAT, LIA2, LIB2, NN, IB1, II 
      REAL(DOUBLE) :: REC, RECS, RECL, WW, B 
      logical rrr
!-----------------------------------------------
      IF (IHSH <= 1) RETURN  
      IF (IA == IB) RETURN  
      IF (IA < IB) THEN 
         IAA = IA 
         IBB = IB 
      ELSE 
         IAA = IB 
         IBB = IA 
      ENDIF 
      IF (.NOT.RECOUPLS0(1,IAA,IBB,IBB,IBB,0)) RETURN  
      rrr = RECOUPLS0(1,IAA,IBB,IBB,IBB,0)
      IF (.NOT.RECOUPLS0(2,IAA,IBB,IBB,IBB,1)) RETURN  
      rrr = RECOUPLS0(2,IAA,IBB,IBB,IBB,1)
      IF (.NOT.RECOUPLS0(3,IAA,IBB,IBB,IBB,1)) RETURN  
      rrr = RECOUPLS0(3,IAA,IBB,IBB,IBB,1)
      CALL RECOUPLS2 (3, IAA, IBB, 1, 0, IAT, REC) 
      IF (IAT == 0) RETURN  
      LIA2 = LJ(IA)*2 
      LIB2 = LJ(IB)*2 
      CALL RECOUPLS2 (2, IAA, IBB, LIB2, 0, IAT, REC) 
      IF (IAT == 0) RETURN  
      CALL RECOUPLS2 (3, IAA, IBB, 1, 1, IAT, RECS) 
      CALL RECOUPLS2 (2, IAA, IBB, LIB2, 1, IAT, RECL) 
      CALL HIBFF (IA, IB, IA, IA, 2) 
      CALL A1A2LS (IK1, IK2, BK1, BK2, ID1, ID2, BD1, BD2, (-HALF), HALF, WW) 
      IF (DABS(WW) > EPS) THEN 
         B = WW*RECL*RECS*HALF*DSQRT(DBLE(4*LJ(IA)+2)) 
         NN = 0 
         IB1 = IBB - 1 
         NN = SUM(NOSH1(IAA:IB1)) 
         IF ((NN/2)*2 == NN) B = -B 
         IF (DABS(B) > EPS) CALL SAVENON (4, B, LJ(IA), 0, IJFUL(IB), 0, IJFUL(&
            IA), JA, JB, 0) 
      ENDIF 
      RETURN  
      END SUBROUTINE LMATRIX2 
