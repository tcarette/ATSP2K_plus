*
*     -------------------------------------------------------------
*      S O O 1 1 2 2  
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
      SUBROUTINE SOO1122(IG,KL,IA,IB,XXX)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/DIAGNL/IDIAG,JA,JB
      EXTERNAL XXX
      L1=LJ(IA)
      L2=LJ(IB)
      CALL SOOA(L1,L1,L2,L2,KL,AP1)
      IF(DABS(AP1).LT.EPS) RETURN
      LA=IJFUL(IA)
      LB=IJFUL(IB)
      KLP=KL+1
      KLM=KL-1
      G1=ZERO
      G2=ZERO
      IF(KLM.GE.0) THEN
        CALL SOOB(1,L1,L2,KLM,KL,A1)
	IF(DABS(A1).GT.EPS) THEN
          CALL XXX(KLM,1,KL,0,1,IA,IB,IB,IB,IB,C1,CC1)
          CALL XXX(KLM,0,KL,1,1,IA,IB,IB,IB,IB,C2,CC2)
          C1=C1*A1
          C2=TWO*C2*A1
          G1=C1+C2
	ENDIF
      ENDIF
      IF(KLP.LE.IG) THEN
        CALL SOOB(1,L1,L2,KLP,KL,A1)
	IF(DABS(A1).GT.EPS) THEN
          CALL XXX(KLP,1,KL,0,1,IA,IB,IB,IB,IB,C3,CC3)
          CALL XXX(KLP,0,KL,1,1,IA,IB,IB,IB,IB,C4,CC4)
          C3=C3*A1
          C4=TWO*C4*A1
          G2=C3+C4
	ENDIF
      ENDIF
      CALL XXX(KL,1,KL,0,1,IA,IB,IB,IB,IB,C55,CC55)
      CALL XXX(KL,0,KL,1,1,IA,IB,IB,IB,IB,C66,CC66)
*
*    ( k k )  tensorine struktura        V   integralas
*
      CALL SOOB(2,L1,L2,KL,KL,A1)
      C7=C55*A1
      C8=TWO*C66*A1
      G4=AP1*(C7+C8)
      IF(DABS(G4).GT.EPS)CALL SAVENON(7,G4,KL-1,LA,LA,LB,LB,JA,JB,0)
*
*    ( k  k+1)  + ( k  k )              
*
      CALL SOOB(1,L1,L2,KL,KL,A1)
      C5=C55*A1
      C6=TWO*C66*A1
      G3=C5+C6
      IF(KL-1.GE.0) THEN
        A1=AP1*(G1+DBLE(KL+1)*G3)
        IF(DABS(A1).GT.EPS)CALL SAVENON(6,A1,KL-2,LA,LA,LB,LB,JA,JB,0)
      ENDIF
*
*    ( k  k-1)  + ( k  k )              
*
      A2=AP1*(G2-DBLE(KL)*G3)
      IF(DABS(A2).GT.EPS)CALL SAVENON(6,A2,KL,LA,LA,LB,LB,JA,JB,0)
      RETURN
      END
