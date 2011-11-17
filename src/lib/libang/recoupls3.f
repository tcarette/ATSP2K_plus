*
*     --------------------------------------------------------------
*     R E C O U P L S 3 
*     --------------------------------------------------------------
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*
      SUBROUTINE RECOUPLS3(K,JA1,JA2,JA3,K1,K2,K3,IRE,IAT,REC)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      IF((JA3.GT.JA1).AND.(JA3.GT.JA2)) THEN
        IF(JA1-JA2.LT.0) THEN
          CALL RECOUPLS31(K,JA1,JA2,JA3,K1,K2,K3,IRE,IAT,REC)
        ELSEIF(JA1-JA2.GT.0) THEN
          CALL RECOUPLS31(K,JA2,JA1,JA3,K2,K1,K3,IRE,IAT,REC)
          IFAZ=K1+K2-K3
          IF((IFAZ/4)*4.NE.IFAZ)REC=-REC
	ELSE
	  STOP
	ENDIF
      ELSEIF((JA3.LT.JA1).AND.(JA3.LT.JA2)) THEN
        IF(JA1-JA2.LT.0) THEN
          CALL RECOUPLS31(K,JA3,JA1,JA2,K3,K1,K2,IRE,IAT,REC)
          IF((K3/2)*2.NE.K3)REC=-REC
        ELSEIF(JA1-JA2.GT.0) THEN
          CALL RECOUPLS31(K,JA3,JA2,JA1,K3,K2,K1,IRE,IAT,REC)
          IFAZ=K1+K2+K3
          IF((IFAZ/4)*4.NE.IFAZ)REC=-REC
	ELSE
	  STOP
	ENDIF
      ELSE
        IF(JA1-JA2.LT.0) THEN
          CALL RECOUPLS31(K,JA1,JA3,JA2,K1,K3,K2,IRE,IAT,REC)
          IFAZ=K1-K2-K3
          IF((IFAZ/4)*4.NE.IFAZ)REC=-REC
        ELSEIF(JA1-JA2.GT.0) THEN
          CALL RECOUPLS31(K,JA2,JA3,JA1,K2,K3,K1,IRE,IAT,REC)
          IF((K1/2)*2.NE.K1)REC=-REC
	ELSE
	  STOP
	ENDIF
      ENDIF
      RETURN
      END
