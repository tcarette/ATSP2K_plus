*
*     -------------------------------------------------------------
*      C 1 E 0 S M
*     -------------------------------------------------------------
*
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
*                                                 ---         ---  *
*                                                 I  Q   1  C   I  *
*     CLEBSCH - GORDAN COEFFICIENT:               I             I  *
*                                                 I  QM  0  CM  I  *
*                                                 ---         ---  *
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vilnius,  Lithuania                          December 1993   *
*
      SUBROUTINE C1E0SM(Q,QM,C,CM,A)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      A=ZERO
      IIQ=TWO*Q+TENTH
      IIC=TWO*C+TENTH
      IF(ITTK(IIQ,IIC,2).EQ.0)RETURN
      IF(ABS(QM-CM).GT.EPS)GO TO 4
      IF((Q+TENTH).LT.ABS(QM))GO TO 4
      IF((C+TENTH).LT.ABS(CM))GO TO 4
      IF(ABS(QM).GT.EPS)GO TO 5
      IS=Q+C+ONE+TENTH
      IF((IS/2)*2.NE.IS)GO TO 4
    5 IG=Q-C+TWO+TENTH
      IF(IG.LE.0)GO TO 4
      IF(IG.GT.3)GO TO 4
      GO TO(2,1,3),IG
    3 A=-SQRT(((C+CM+ONE)*(C-CM+ONE))/((C+ONE)*(TWO*C+THREE)))
      RETURN
    1 A=CM/SQRT(C*(C+ONE))
      RETURN
    2 A=SQRT(((C+CM)*(C-CM))/((TWO*C-ONE)*C))
      RETURN
    4 A=ZERO
      RETURN
      END
