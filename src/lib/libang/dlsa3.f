*
*     --------------------------------------------------------------
*      D L S A 3
*     --------------------------------------------------------------
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*
      SUBROUTINE DLSA3(K,JA1,JA2,KA,IRE,IAT,REC)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      REC=ZERO
      AA=ONE
      IMAX=JA2-1
      IMIN=JA1+1
      IF(JA1.EQ.1)IMIN=IMIN+1
      IF(IMIN.LT.JA2) THEN
        DO 1 I=IMIN,IMAX
          JI=J1QN1(I,K)-1
          K1=IHSH+I-2
          ITI=J1QN1(K1,K)-1
          ITIS=J1QN2(K1,K)-1
          K1=K1+1
          ITI1=J1QN1(K1,K)-1
          ITI1S=J1QN2(K1,K)-1
          IF(IRE.EQ.0) THEN
            IF(IXJTIK(KA,ITIS,ITI,JI,ITI1,ITI1S).EQ.0)RETURN
          ELSE
            CALL SIXJ(KA,ITIS,ITI,JI,ITI1,ITI1S,0,A3)
            A3=A3*DSQRT(DBLE((ITI+1)*(ITI1S+1)))
            IFAZ=KA+JI+ITI+ITI1S
            IF((IFAZ/4)*4.NE.IFAZ)A3=-A3
            AA=AA*A3
          ENDIF
    1   CONTINUE
      ENDIF
      REC=AA
      IAT=1
      RETURN
      END
