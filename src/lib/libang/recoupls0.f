*
*     --------------------------------------------------------------
*     R E C O U P L S 0
*     --------------------------------------------------------------
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*
      LOGICAL FUNCTION RECOUPLS0(K,JA1,JA2,JA3,JA4,KA)
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      RECOUPLS0=.TRUE.
      IF(IHSH.EQ.1)RETURN
      IF(JA1.EQ.1.AND.JA2.EQ.2) GO TO 1
      IF(KA.NE.0)GO TO 5
*
*  CASES WHEN :          KA = 0
*                  OR    JA1 = JA2
*                  OR    JA1 = 1    JA2 = 2
*
    1 DO 3 I=1,IHSH
      IJ=IHSH+I-1
      IF(I.EQ.1)GO TO 4
      IF(J1QN1(IJ,K).NE.J1QN2(IJ,K)) RECOUPLS0=.FALSE.
    4 IF(KA.EQ.0)GO TO 9
      IF(I.EQ.JA1)GO TO 3
      IF(I.EQ.JA2)GO TO 3
    9 CONTINUE
      IF(I.EQ.JA1.AND.K.EQ.1)GO TO 3
      IF(I.EQ.JA2.AND.K.EQ.1)GO TO 3
      IF(I.EQ.JA3.AND.K.EQ.1)GO TO 3
      IF(I.EQ.JA4.AND.K.EQ.1)GO TO 3
      IF(J1QN1(I,K).NE.J1QN2(I,K)) RECOUPLS0=.FALSE.
    3 CONTINUE
      RETURN
*
*  OTHER CASES
*
    5 CONTINUE
      IA1=JA1-1
      IA2=JA2-1
      IF(JA1.EQ.1)IA1=JA1
      DO 6 I=1,IHSH
      IJ=IHSH+I
      IF(I.EQ.IHSH)GO TO 7
      IF(I.GE.IA1.AND.I.LT.IA2)GO TO 7
      IF(J1QN1(IJ,K).NE.J1QN2(IJ,K)) RECOUPLS0=.FALSE.
    7 IF(I.EQ.JA1)GO TO 6
      IF(I.EQ.JA2)GO TO 6
      IF((KA.EQ.2).AND.(I.EQ.JA3))GO TO 6
      IF((KA.EQ.3).AND.(I.EQ.JA3))GO TO 6
      IF((KA.EQ.3).AND.(I.EQ.JA4))GO TO 6
      IF(J1QN1(I,K).NE.J1QN2(I,K)) RECOUPLS0=.FALSE.
    6 CONTINUE
      RETURN
      END
