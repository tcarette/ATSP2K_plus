*
*     --------------------------------------------------------------
*     R L S P 3 1
*     --------------------------------------------------------------
*                                                                  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE RLSP31(K,JA1,JA2,JA3,K1,K2,K3,K4,KA,IRE,IAT,REC)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
*
      REC=ONE
      IAT=1
      IF((IHSH-JA3).GT.1) THEN
        IAT=0
        CALL DLSA3(K,JA3,IHSH,KA,IRE,IAT,RE)
        IF(IAT.EQ.0)RETURN
        REC=RE*REC
      ENDIF
*
      IF(JA3.NE.IHSH) THEN
        IAT=0
        CALL DLSA5(K,JA3,KA,IRE,IAT,RE)
        IF(IAT.EQ.0)RETURN
        REC=RE*REC
      ENDIF
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
