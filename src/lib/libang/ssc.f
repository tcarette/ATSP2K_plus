*
*     -------------------------------------------------------------
*      S S C
*     -------------------------------------------------------------
*                                                                  *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
*                                                                  *
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vanderbilt University Nashville, USA           October 1995  *
*
      SUBROUTINE SSC(IG,KL,IA,IB,IC,ID,IIA,IIB,IIC,IID,IREZ,XXX)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/DIAGNL/IDIAG,JA,JB
      EXTERNAL XXX
      KLP=KL+2
C      IF(IG-1.LT.KLP) RETURN
      L1=LJ(IIA)
      L2=LJ(IIB)
      L3=LJ(IIC)
      L4=LJ(IID)
      CALL SSA(L1,L2,L3,L4,KLP,KL,A1)
      CALL SSA(L2,L1,L4,L3,KLP,KL,A2)
      IF((DABS(A1)+DABS(A2)).LT.TWO*EPS) RETURN
      LA=IJFUL(IIA)
      LB=IJFUL(IIB)
      LC=IJFUL(IIC)
      LD=IJFUL(IID)
      CALL XXX(KLP,1,KL,1,2,IA,IB,IC,ID,IREZ,C,CC)
      C=C*A1
      CC=CC*A2
      IF(DABS(C).GT.EPS)CALL SAVENON(8,C,KL,LA,LB,LC,LD,JA,JB,0)
      IF(DABS(CC).GT.EPS)CALL SAVENON(8,CC,KL,LB,LA,LD,LC,JA,JB,0)
      RETURN
      END
