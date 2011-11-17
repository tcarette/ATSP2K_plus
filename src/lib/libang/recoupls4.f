*
*     --------------------------------------------------------------
*     R E C O U P L S 4
*     --------------------------------------------------------------
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*
      SUBROUTINE RECOUPLS4(K,JA1,JA2,JA3,JA4,K1,K2,K3,K4,KA,
     *IRE,IAT,REC)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      S1=DBLE(J1QN1(JA1,K))
      S2=DBLE(J1QN1(JA2,K))
      S3=DBLE(J1QN1(JA3,K))
      S4=DBLE(J1QN1(JA4,K))
      S=S1*S2*S3*S4
      REC=ONE/DSQRT(S)
      IA4=J1QN1(JA4,K)-1
      IB4=J1QN2(JA4,K)-1
      REC=REC*DSQRT(DBLE(IA4+1))/DSQRT(DBLE((K4+1)*(IB4+1)))
*
      ISKR=JA3-JA2
      IF(ISKR.GT.1) THEN
        IAT=0
        CALL DLSA3(K,JA2,JA3,KA,IRE,IAT,RE)
        IF(IAT.EQ.0)RETURN
        REC=RE*REC
      ENDIF
*
      ISKR=JA4-JA3
      IF(ISKR.GT.1) THEN
        IAT=0
        CALL DLSA3(K,JA3,JA4,K4,IRE,IAT,RE)
        IF(IAT.EQ.0)RETURN
        REC=RE*REC
      ENDIF
*
      IAT=0
      CALL DLSA2(K,JA1,JA4,K4,IRE,IAT,RE)
      IF(IAT.EQ.0)RETURN
      REC=RE*REC
*
      IAT=0
      CALL DLSA4(K,JA1,JA2,K1,K2,KA,IRE,IAT,RE)
      IF(IAT.EQ.0)RETURN
      REC=RE*REC
*
      IAT=0
      CALL DLSA4(K,JA2,JA3,KA,K3,K4,IRE,IAT,RE)
      IF(IAT.EQ.0)RETURN
      REC=RE*REC
      IF(JA1.EQ.1.AND.JA2.EQ.2)RETURN
*
      IAT=0
      CALL DLSA1(K,JA1,K1,IRE,IAT,RE)
      IF(IAT.EQ.0)RETURN
      REC=RE*REC
*
      ISKR=JA2-JA1
      IF(JA1.EQ.1)ISKR=JA2-1-JA1
      IF(ISKR.LE.1)RETURN
      IAT=0
      CALL DLSA3(K,JA1,JA2,K1,IRE,IAT,RE)
      REC=RE*REC
      RETURN
      END
