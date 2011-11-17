*
*     --------------------------------------------------------------
*     R E C O U P L S 2 
*     --------------------------------------------------------------
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*
      SUBROUTINE RECOUPLS2(K,JA1,JA2,KA,IRE,IAT,REC)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      IAT=0
      S=DBLE(J1QN1(JA1,K))
      SS=DBLE(J1QN1(JA2,K))
      REC=ONE/DSQRT(S*SS)
      IF((IRE.NE.0).AND.(KA.EQ.0)) THEN
        IAT=1
      ELSE
        IA1=J1QN1(JA1,K)-1
        IB1=J1QN2(JA1,K)-1
        IA2=J1QN1(JA2,K)-1
        IB2=J1QN2(JA2,K)-1
        IAT=0
        CALL DLSA2(K,JA1,JA2,KA,IRE,IAT,RE)
        IF(IAT.EQ.0)RETURN
        REC=RE*REC*DSQRT(DBLE(IA2+1))/DSQRT(DBLE((KA+1)*(IB2+1)))
        IF(JA1.EQ.1.AND.JA2.EQ.2)RETURN
        IAT=0
        CALL DLSA1(K,JA1,KA,IRE,IAT,RE)
        IF(IAT.EQ.0)RETURN
        REC=RE*REC
        ISKR=JA2-JA1
        IF(JA1.EQ.1)ISKR=JA2-1-JA1
        IF(ISKR.LE.1)RETURN
        IAT=0
        CALL DLSA3(K,JA1,JA2,KA,IRE,IAT,RE)
        REC=RE*REC
      ENDIF
      RETURN
      END
