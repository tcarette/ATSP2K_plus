*
*     -------------------------------------------------------------
*      S O O 1 1 1 1 
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
      SUBROUTINE SOO1111(IG,KL,IA)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/DIAGNL/IDIAG,JA,JB
      KLP=KL+1
      KLM=KL-1
      LA=IJFUL(IA)
      IF(KLM.GE.0) THEN
        CALL SOO1(LJ(IA),KLM,1,KL,0,A1)
        IF(DABS(A1).GT.EPS) THEN
          CALL TWO11(KLM,1,KL,0,1,IA,C1)
	  C1=C1*A1
        ELSE
          C1=ZERO
        ENDIF
        CALL SOO1(LJ(IA),KLM,0,KL,1,A1)
        IF(DABS(A1).GT.EPS) THEN
          CALL TWO11(KLM,0,KL,1,1,IA,C2)
	  C2=C2*A1
        ELSE
          C2=ZERO
        ENDIF
        C=C1+C2
        IF(DABS(C).GT.EPS)CALL SAVENON(6,C,KLM-1,LA,LA,LA,LA,JA,JB,0)
      ENDIF
      IF(KLP.LE.IG) THEN
        CALL SOO1(LJ(IA),KLP,1,KL,0,A1)
        IF(DABS(A1).GT.EPS) THEN
          CALL TWO11(KLP,1,KL,0,1,IA,C1)
	  C1=C1*A1
        ELSE
          C1=ZERO
        ENDIF
        CALL SOO1(LJ(IA),KLP,0,KL,1,A1)
        IF(DABS(A1).GT.EPS) THEN
          CALL TWO11(KLP,0,KL,1,1,IA,C2)
	  C2=C2*A1
        ELSE
          C2=ZERO
        ENDIF
        C=C1+C2
        IF(DABS(C).GT.EPS)CALL SAVENON(6,C,KL,LA,LA,LA,LA,JA,JB,0)
      ENDIF
      RETURN
      END
