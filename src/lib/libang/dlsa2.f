*
*     --------------------------------------------------------------
*      D L S A 2
*     --------------------------------------------------------------
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*
      SUBROUTINE DLSA2(K,JA1,JA2,KA,IRE,IAT,REC)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      REC=ZERO
      IA1=J1QN1(JA1,K)-1
      IB1=J1QN2(JA1,K)-1
      IA2=J1QN1(JA2,K)-1
      IB2=J1QN2(JA2,K)-1
      IF(JA1.EQ.1.AND.JA2.EQ.2) THEN
         IT2=IA1
         IT2S=IB1
         N1=IHSH+1
         J2=J1QN1(N1,K)-1
      ELSE
         N1=IHSH+JA2-1
         J2=J1QN1(N1,K)-1
         N2=IHSH+JA2-2
         IT2=J1QN1(N2,K)-1
         IT2S=J1QN2(N2,K)-1
      ENDIF
      IF(IRE.EQ.0) THEN
         IF(IXJTIK(KA,IB2,IA2,J2,IT2,IT2S).EQ.0)RETURN
         IAT=1
      ELSE
         CALL SIXJ(KA,IB2,IA2,J2,IT2,IT2S,0,A2)
         REC=A2*DSQRT(DBLE((IB2+1)*(IT2+1)))
         IFAZ=J2+IT2S+IA2+KA
         IF((IFAZ/4)*4.NE.IFAZ)REC=-REC
         IAT=1
      ENDIF
      RETURN
      END
