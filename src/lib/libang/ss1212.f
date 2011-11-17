*
*     -------------------------------------------------------------
*      S S 1 2 1 2 
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
C      CALL COULOMBLS(LJ(IA),LJ(IB),LJ(IA),LJ(IB),KL1,A1)
C      A1=A1*TWO*DSQRT(DBLE(2*KL1+1))
CGG elektrostatine
      SUBROUTINE SS1212(IG,KL,IA,IB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/DIAGNL/IDIAG,JA,JB
      KLP=KL+2
C      IF(IG-1.LT.KLP) RETURN
      LA=IJFUL(IA)
      LB=IJFUL(IB)
      CALL SS1(LJ(IA),LJ(IB),KLP,KL,A1)
      IF(DABS(A1).GT.EPS) THEN
        CALL TWO12(KLP,1,KL,1,2,IA,IB,C)
	C=C*A1
        IF(DABS(C).GT.EPS)CALL SAVENON(8,C,KL,LA,LB,LA,LB,JA,JB,0)
      ENDIF
      CALL SS1(LJ(IB),LJ(IA),KLP,KL,A1)
      IF(DABS(A1).LE.EPS) RETURN
      CALL TWO12(KL,1,KLP,1,2,IA,IB,C)
      C=C*A1
      IF(DABS(C).GT.EPS)CALL SAVENON(8,C,KL,LB,LA,LB,LA,JA,JB,0)
      RETURN
      END
