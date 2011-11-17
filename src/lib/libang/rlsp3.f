*
*
*     --------------------------------------------------------------
*     R L S P 3 
*     --------------------------------------------------------------
*                                                                  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE RLSP3(K,JA1,JA2,JA3,K1,K2,K3,K4,KA,REC)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      REC=ZERO
*
      IF((JA3.GT.JA1).AND.(JA3.GT.JA2)) THEN
        IF(JA1-JA2.LT.0) THEN
* 1 2 3
          CALL RLSP31(K,JA1,JA2,JA3,K1,K2,K3,K4,KA,0,IAT,R)
          IF(IAT.EQ.0) RETURN
          CALL RLSP32(K,JA1,JA2,JA3,K1,K2,K3,K4,KA,0,IAT,R)
          IF(IAT.EQ.0) RETURN
          CALL RLSP31(K,JA1,JA2,JA3,K1,K2,K3,K4,KA,1,IAT,REC)
          CALL RLSP32(K,JA1,JA2,JA3,K1,K2,K3,K4,KA,1,IAT,R1)
	  REC=REC*R1
        ELSEIF(JA1-JA2.GT.0) THEN
* 2 1 3
          CALL RLSP31(K,JA2,JA1,JA3,K2,K1,K3,K4,KA,0,IAT,R)
          IF(IAT.EQ.0) RETURN
          CALL RLSP32(K,JA2,JA1,JA3,K2,K1,K3,K4,KA,0,IAT,R)
          IF(IAT.EQ.0) RETURN
          CALL RLSP31(K,JA2,JA1,JA3,K2,K1,K3,K4,KA,1,IAT,REC)
          CALL RLSP32(K,JA2,JA1,JA3,K2,K1,K3,K4,KA,1,IAT,R1)
	  REC=REC*R1
          IF(MOD(K1+K2-K3,4).NE.0)REC=-REC
	ELSE
	  STOP
	ENDIF
      ELSEIF((JA3.LT.JA1).AND.(JA3.LT.JA2)) THEN
        IF(JA1-JA2.LT.0) THEN
* 3 1 2
          IPL=ITREXG(K1,K4,K2,KA,IKKL)+1
          IF(IKKL.EQ.0) RETURN
          IGL=IPL+IKKL-1
          CALL RLSP31(K,JA3,JA1,JA2,K4,K1,J12,K2,KA,0,IAT,R)
          IF(IAT.EQ.0) RETURN
          DO 3 I=IPL,IGL,2
            J12=I-1
            CALL DLSA6(K,K4,K3,KA,K1,K2,J12,0,IAT,R)
            IF(IAT.NE.0) THEN
              CALL RLSP32(K,JA3,JA1,JA2,K4,K1,J12,K2,KA,0,IAT,R)
              IF(IAT.NE.0) THEN
                CALL RLSP32(K,JA3,JA1,JA2,K4,K1,J12,K2,KA,1,IAT,R1)
                CALL DLSA6(K,K4,K3,KA,K1,K2,J12,1,IAT,R2)
                IF(MOD(2*K1+K2+K4-J12-K3,4).NE.0)R1=-R1
                REC=REC+R1*R2
              ENDIF
            ENDIF
    3     CONTINUE
	  IF(DABS(REC).LT.EPS)RETURN
          CALL RLSP31(K,JA3,JA1,JA2,K4,K1,J12,K2,KA,1,IAT,R1)
          REC=REC*R1
        ELSEIF(JA1-JA2.GT.0) THEN
* 3 2 1
          IPL=ITREXG(K2,K4,K1,KA,IKKL)+1
          IF(IKKL.EQ.0) RETURN
          IGL=IPL+IKKL-1
          CALL RLSP31(K,JA3,JA2,JA1,K4,K2,J12,K1,KA,0,IAT,R)
          IF(IAT.EQ.0) RETURN
          DO 4 I=IPL,IGL,2
            J12=I-1
            CALL DLSA6(K,K4,K3,KA,K2,K1,J12,0,IAT,R)
            IF(IAT.NE.0) THEN
              CALL RLSP32(K,JA3,JA2,JA1,K4,K2,J12,K1,KA,0,IAT,R)
              IF(IAT.NE.0) THEN
                CALL RLSP32(K,JA3,JA2,JA1,K4,K2,J12,K1,KA,1,IAT,R1)
                CALL DLSA6(K,K4,K3,KA,K2,K1,J12,1,IAT,R2)
                IF(MOD(3*K2+2*K1+K4-J12-2*K3,4).NE.0)R1=-R1
                REC=REC+R1*R2
              ENDIF
            ENDIF
    4     CONTINUE
	  IF(DABS(REC).LT.EPS)RETURN
          CALL RLSP31(K,JA3,JA2,JA1,K4,K2,J12,K1,KA,1,IAT,R1)
          REC=REC*R1
	ELSE
	  STOP
	ENDIF
      ELSE
        IF(JA1-JA2.LT.0) THEN
* 1 3 2
          IPL=ITREXG(K1,K4,K2,KA,IKKL)+1
          IF(IKKL.EQ.0) RETURN
          IGL=IPL+IKKL-1
          CALL RLSP31(K,JA1,JA3,JA2,K1,K4,J12,K2,KA,0,IAT,R)
          IF(IAT.EQ.0) RETURN
          DO 5 I=IPL,IGL,2
            J12=I-1
            CALL DLSA6(K,K4,K3,KA,K1,K2,J12,0,IAT,R)
            IF(IAT.NE.0) THEN
              CALL RLSP32(K,JA1,JA3,JA2,K1,K4,J12,K2,KA,0,IAT,R)
              IF(IAT.NE.0) THEN
                CALL RLSP32(K,JA1,JA3,JA2,K1,K4,J12,K2,KA,1,IAT,R1)
                CALL DLSA6(K,K4,K3,KA,K1,K2,J12,1,IAT,R2)
                REC=REC+R1*R2
              ENDIF
            ENDIF
    5     CONTINUE
	  IF(DABS(REC).LT.EPS)RETURN
          CALL RLSP31(K,JA1,JA3,JA2,K1,K4,J12,K2,KA,1,IAT,R1)
          REC=REC*R1
          IF(MOD(K1+K2-K3,4).NE.0)REC=-REC
        ELSEIF(JA1-JA2.GT.0) THEN
* 2 3 1
          IPL=ITREXG(K2,K4,K1,KA,IKKL)+1
          IF(IKKL.EQ.0) RETURN
          IGL=IPL+IKKL-1
          CALL RLSP31(K,JA2,JA3,JA1,K2,K4,J12,K1,KA,0,IAT,R)
          IF(IAT.EQ.0) RETURN
          DO 6 I=IPL,IGL,2
            J12=I-1
            CALL DLSA6(K,K4,K3,KA,K2,K1,J12,0,IAT,R)
            IF(IAT.NE.0) THEN
              CALL RLSP32(K,JA2,JA3,JA1,K2,K4,J12,K1,KA,0,IAT,R)
              IF(IAT.NE.0) THEN
                CALL RLSP32(K,JA2,JA3,JA1,K2,K4,J12,K1,KA,1,IAT,R1)
                CALL DLSA6(K,K4,K3,KA,K2,K1,J12,1,IAT,R2)
                REC=REC+R1*R2
              ENDIF
            ENDIF
    6     CONTINUE
	  IF(DABS(REC).LT.EPS)RETURN
          CALL RLSP31(K,JA2,JA3,JA1,K2,K4,J12,K1,KA,1,IAT,R1)
          REC=REC*R1
	ELSE
	  STOP
	ENDIF
      ENDIF
      REC=REC/DSQRT(DBLE(J1QN1(JA1,K)*J1QN1(JA2,K)*J1QN1(JA3,K)))
      RETURN
      END
