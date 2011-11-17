*
*     -------------------------------------------------------------
*      R L S P 0
*     -------------------------------------------------------------
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE RLSP0(K,JA1,JA2,KA,IAT)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      IF(K.NE.1) THEN
        IAT=0
        KK=IHSH+IHSH-1
        LK=J1QN1(KK,K)-1
        LD=J1QN2(KK,K)-1
        IF(ITTK(LK,LD,KA).EQ.0)RETURN
      ENDIF
      IAT=1
      IF(IHSH.EQ.1)RETURN
      DO 1 I=1,IHSH
        IF(JA1.NE.I) THEN
          IF(JA2.NE.I) THEN
            IF(J1QN1(I,K).NE.J1QN2(I,K))IAT=0
          ENDIF
        ENDIF
    1 CONTINUE
      IF(IAT.EQ.0)RETURN
      IF(K.EQ.1)RETURN
      IF(IHSH.LE.2)RETURN
      IF(JA1.LE.2)RETURN
      DO 2 J=3,JA1
        JJ=IHSH-2+J
        IF(J1QN1(JJ,K).NE.J1QN2(JJ,K))IAT=0
        IF(IAT.EQ.0)RETURN
    2 CONTINUE
C cia dabar
      ISKR=IHSH-JA2
      IF(ISKR.GT.0) THEN
	DO 3 JI=1,ISKR
	  KK=IHSH+JA2-2+JI
          LK=J1QN1(KK,K)-1
          LD=J1QN2(KK,K)-1
          IF(ITTK(LK,LD,KA).EQ.0)IAT=0
          IF(IAT.EQ.0)RETURN
    3   CONTINUE
      ENDIF
C cia dabar
      RETURN
      END
