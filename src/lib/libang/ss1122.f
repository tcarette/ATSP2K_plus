*
*     -------------------------------------------------------------
*      S S 1 1 2 2
*     -------------------------------------------------------------
*                                                                  *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
*                                                                  *
*                                                                  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
*
CGG elektrostatine
C      CALL COULOMBLS(LJ(IA),LJ(IA),LJ(IB),LJ(IB),KL1,A1)
C      A1=A1*TWO*DSQRT(DBLE(2*KL1+1))
CGG elektrostatine
      SUBROUTINE SS1122(IG,KL,IA,IB)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/DIAGNL/IDIAG,JA,JB
      KLP=KL+2
C      IF(IG-1.LT.KLP) RETURN
      LA=IJFUL(IA)
      LB=IJFUL(IB)
      CALL SSA(LJ(IA),LJ(IA),LJ(IB),LJ(IB),KLP,KL,A1) 
      IF(DABS(A1).LT.EPS) RETURN
      INE3=0
      CALL TWO2(KLP,1,KL,1,2,IA,IB,IB,IB,INE3,C,CC)
      C=C*A1
      IF(DABS(C).GT.EPS)CALL SAVENON(8,C,KL,LA,LA,LB,LB,JA,JB,0)
      RETURN
      END
