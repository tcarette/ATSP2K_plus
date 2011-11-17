*
*     -------------------------------------------------------------
*      S O O 1 2 1 2
*     -------------------------------------------------------------
*                                                                  *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
*
      SUBROUTINE SOO1212(IG,KL,IA,IB)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/DIAGNL/IDIAG,JA,JB
      L1=LJ(IA)
      L2=LJ(IB)
      CALL SOOA(L1,L2,L1,L2,KL,AP1)
      IF(DABS(AP1).LT.EPS) RETURN
      LB=IJFUL(IB)
      LA=IJFUL(IA)
      KLP =KL+1
      KLM =KL-1
      IF(KLM.GE.0) THEN
      CALL SOOB(1,L1,L1,KLM,KL,A1)
        IF(DABS(A1).GT.EPS) THEN
          CALL TWO12(KLM,1,KL,0,1,IA,IB,C1)
          CALL TWO12(KLM,0,KL,1,1,IA,IB,C2)
          C=C1+TWO*C2
          C=AP1*A1*C
          IF(DABS(C).GT.EPS)CALL SAVENON(6,C,KLM-1,LB,LA,LB,LA,JA,JB,0)
        ENDIF
*
        CALL SOOB(1,L2,L2,KLM,KL,A1)
        IF(DABS(A1).GT.EPS) THEN
          CALL TWO12(KL,0,KLM,1,1,IA,IB,C1)
          CALL TWO12(KL,1,KLM,0,1,IA,IB,C2)
          C=C1+TWO*C2
          C=AP1*A1*C
          IF(DABS(C).GT.EPS)CALL SAVENON(6,C,KLM-1,LA,LB,LA,LB,JA,JB,0)
        ENDIF
      ENDIF
      IF(KLP.LE.IG) THEN
        CALL SOOB(1,L1,L1,KLP,KL,A1)
        IF(DABS(A1).GT.EPS) THEN
          CALL TWO12(KLP,1,KL,0,1,IA,IB,C3)
          CALL TWO12(KLP,0,KL,1,1,IA,IB,C4)
          C=C3+TWO*C4
          C=AP1*A1*C
          IF(DABS(C).GT.EPS)CALL SAVENON(6,C,KL,LA,LB,LA,LB,JA,JB,0)
        ENDIF
*
        CALL SOOB(1,L2,L2,KLP,KL,A1)
	IF(DABS(A1).GT.EPS) THEN
          CALL TWO12(KL,0,KLP,1,1,IA,IB,C3)
          CALL TWO12(KL,1,KLP,0,1,IA,IB,C4)
          C=C3+TWO*C4
          C=AP1*A1*C
          IF(DABS(C).GT.EPS)CALL SAVENON(6,C,KL,LB,LA,LB,LA,JA,JB,0)
	ENDIF
      ENDIF
      RETURN
      END
