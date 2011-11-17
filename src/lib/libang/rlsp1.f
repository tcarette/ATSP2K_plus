*
*     -------------------------------------------------------------
*      R L S P 1
*     -------------------------------------------------------------
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE RLSP1(K,JA1,KA,IRE,IAT,REC)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      IAT=1
      REC=ONE/SQRT(DBLE(J1QN1(JA1,K)))
      IF(IHSH.EQ.1)RETURN
      IAT=0
      IF(IRE.NE.0) THEN
        IF(KA.EQ.0) THEN
          IAT=1
          RETURN
        ENDIF
      ENDIF
      IF(IHSH.NE.2) THEN
        CALL DLSA5(K,JA1,KA,IRE,IAT,RE)
        REC=RE*REC
        IF(IAT.EQ.0)RETURN
        IF(JA1.EQ.IHSH)RETURN
        IAT=0
      ENDIF
      CALL DLSA1(K,JA1,KA,IRE,IAT,RE)
      IF(IAT.EQ.0)RETURN
      REC=RE*REC
      IF(IHSH.EQ.2)RETURN
      ISKR=IHSH-JA1
      IF(JA1.EQ.1)ISKR=IHSH-1-JA1
      IF(ISKR.LE.1)RETURN
      IAT=0
      CALL DLSA3(K,JA1,IHSH,KA,IRE,IAT,RE)
      REC=RE*REC
      RETURN
      END
