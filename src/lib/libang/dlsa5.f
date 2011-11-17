*
*     -------------------------------------------------------------
*      D L S A 5
*     -------------------------------------------------------------
*                                                                  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE DLSA5(K,JA1,KA,IRE,IAT,REC)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      REC=ZERO
      K1=IHSH+IHSH-1
      ITI1=J1QN1(K1,K)-1
      ITI1S=J1QN2(K1,K)-1
      K1=K1-1
      IF(JA1.EQ.IHSH) THEN
        ITI=J1QN1(IHSH,K)-1
        ITIS=J1QN2(IHSH,K)-1
        JI=J1QN1(K1,K)-1
      ELSE
        JI=J1QN1(IHSH,K)-1
        ITI=J1QN1(K1,K)-1
        ITIS=J1QN2(K1,K)-1
      ENDIF
      IF(IRE.EQ.0) THEN
        IF(IXJTIK(KA,ITIS,ITI,JI,ITI1,ITI1S).NE.0) IAT=1
      ELSE
        CALL SIXJ(KA,ITIS,ITI,JI,ITI1,ITI1S,0,A3)
        REC=A3*DSQRT(DBLE((ITI+1)*(ITI1S+1)))
        IF(MOD(KA+JI+ITIS+ITI1,4).NE.0)REC=-REC
        IAT=1
        IF(JA1.EQ.IHSH)RETURN
        IF(MOD(ITI+ITIS-ITI1S-ITI1+2*JI,4).NE.0)REC=-REC
      ENDIF
      RETURN
      END
