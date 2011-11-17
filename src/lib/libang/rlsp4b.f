*
*     --------------------------------------------------------------
*     R L S P 4 B
*     --------------------------------------------------------------
*                                                                  *
*
      SUBROUTINE RLSP4B(K,JA3,JA4,K5,K6,KA,IRE,IAT,REC)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      REC=ONE/DSQRT(DBLE(J1QN1(JA4,K)))
*
      ISKR=JA4-JA3
      IF(ISKR.GT.1) THEN
        IAT=0
        CALL DLSA3(K,JA3,JA4,K5,IRE,IAT,RE)
        IF(IAT.EQ.0)RETURN
        REC=RE*REC
      ENDIF
*
      ISKR=IHSH-JA4
      IF(ISKR.GT.1) THEN
        IAT=0
        CALL DLSA3(K,JA4,IHSH,KA,IRE,IAT,RE)
        IF(IAT.EQ.0)RETURN
        REC=RE*REC
      ENDIF
*
      IF(JA4.NE.IHSH) THEN
        IAT=0
        CALL DLSA5(K,JA4,KA,IRE,IAT,RE)
        IF(IAT.EQ.0)RETURN
        REC=RE*REC
      ENDIF
      IF(JA3.EQ.JA4) RETURN
*
      IAT=0
      CALL DLSA4(K,JA3,JA4,K5,K6,KA,IRE,IAT,RE)
      IF(IAT.EQ.0)RETURN
      REC=RE*REC
      RETURN
      END
