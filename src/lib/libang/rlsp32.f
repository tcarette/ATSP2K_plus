*
*     --------------------------------------------------------------
*     R L S P 3 2
*     --------------------------------------------------------------
*                                                                  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE RLSP32(K,JA1,JA2,JA3,K1,K2,K3,K4,KA,IRE,IAT,REC)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
*
      REC=ONE
      IF((JA3-JA2).GT.1) THEN
        IAT=0
        CALL DLSA3(K,JA2,JA3,K3,IRE,IAT,RE)
        IF(IAT.EQ.0)RETURN
        REC=RE*REC
      ENDIF
*
      IAT=0
      CALL DLSA4(K,JA1,JA2,K1,K2,K3,IRE,IAT,RE)
      IF(IAT.EQ.0)RETURN
      REC=RE*REC
*
      IAT=0
      CALL DLSA4(K,JA2,JA3,K3,K4,KA,IRE,IAT,RE)
      IF(IAT.EQ.0)RETURN
      REC=RE*REC
      RETURN
      END
