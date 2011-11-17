*
*     -------------------------------------------------------------
*      S O O 1 1 1 1 P
*     -------------------------------------------------------------
*                                                                  *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
*                                                                  *
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Universite Libre de Bruxelles, Belgium         October 1995  *
*
      SUBROUTINE SOO1111P(IG,KL,IA)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/DIAGNL/IDIAG,JA,JB
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
        G1=C1+C2
      ELSE
        G1=ZERO
      ENDIF
      KLP=KL+1
      IF(KLP.LE.IG) THEN
        CALL SOO1(LJ(IA),KLP,1,KL,0,A1)
        IF(DABS(A1).GT.EPS) THEN
          CALL TWO11(KLP,1,KL,0,1,IA,C3)
	  C3=C3*A1
        ELSE
	  C3=ZERO
        ENDIF
        CALL SOO1(LJ(IA),KLP,0,KL,1,A1)
        IF(DABS(A1).GT.EPS) THEN
          CALL TWO11(KLP,0,KL,1,1,IA,C4)
	  C4=C4*A1
        ELSE
	  C4=ZERO
        ENDIF
        G2=C3+C4
      ELSE
        G2=ZERO
      ENDIF
      IF(KL.GT.0) THEN
        CALL SOO1(LJ(IA),KL,1,KL,0,A1)
        IF(DABS(A1).GT.EPS) THEN
          CALL TWO11(KL,1,KL,0,1,IA,C5)
          C5=C5*A1
        ELSE
          C5=ZERO
        ENDIF
        CALL SOO1(LJ(IA),KL,0,KL,1,A1)
        IF(DABS(A1).GT.EPS) THEN
          CALL TWO11(KL,0,KL,1,1,IA,C6)
          C6=C6*A1
        ELSE
          C6=ZERO
        ENDIF
      ELSE
       C5=ZERO
       C6=ZERO
      ENDIF
*
*     ( k k )   tenzorine struktura    V integralas
*
      G3=C5+C6
      G4=TWO*G3
      IF(DABS(G4).GT.EPS)CALL SAVENON(7,G4,KL-1,LA,LA,LA,LA,JA,JB,0)
*
*     ( k k+1 )  +  ( k k ) 
*
      IF(KL-1.GE.0) THEN
        A1=G1+DBLE(KL+1)*G3
        IF(DABS(A1).GT.EPS)CALL SAVENON(6,A1,KL-2,LA,LA,LA,LA,JA,JB,0)
      ENDIF
*
*     ( k k-1 )  +  ( k k ) 
*
      IF(KLP.LE.IG) THEN
        A1=G2-DBLE(KL)*G3
        IF(DABS(A1).GT.EPS)CALL SAVENON(6,A1,KL,LA,LA,LA,LA,JA,JB,0)
      ENDIF
      RETURN
      END
