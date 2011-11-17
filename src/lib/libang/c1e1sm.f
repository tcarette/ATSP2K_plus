*
*     -------------------------------------------------------------
*      C 1 E 1 S M
*     -------------------------------------------------------------
*
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
*                                                 ---         ---  *
*                                                 I  Q   1  C   I  *
*     CLEBSCH - GORDAN COEFFICIENT:               I             I  *
*                                                 I  QM  1  CM  I  *
*                                                 ---         ---  *
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vilnius,  Lithuania                          December 1993   *
*
      SUBROUTINE C1E1SM(Q,QM,SM,C,CM,A)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DIMENSION GC(2)
      GC(1)=ONE
      GC(2)=-ONE
      A=ZERO
      IIQ=TWO*Q+TENTH
      IIC=TWO*C+TENTH
      IF(ITTK(IIQ,IIC,2).EQ.0)RETURN
      IF(ABS(QM+SM-CM).GT.EPS)RETURN
      IF((Q+TENTH).LT.ABS(QM))RETURN
      IF((C+TENTH).LT.ABS(CM))RETURN
      IE=0
      IF(ABS(SM-ONE).LT.EPS)IE=1
      IF(ABS(SM+ONE).LT.EPS)IE=2
      IF(IE.EQ.0)RETURN
      IF(ABS(Q+ONE-C).LT.EPS)GO TO 1
      IF(ABS(Q-C).LT.EPS)GO TO 2
      IF(ABS(Q-ONE-C).GT.EPS)RETURN
      A=SQRT((C-GC(IE)*CM+ONE)*(C-GC(IE)*CM+TWO)/
     *((TWO*C+TWO)*(TWO*C+THREE)))
      RETURN
    1 A=SQRT((C+GC(IE)*CM-ONE)*(C+GC(IE)*CM)/
     *((TWO*C-ONE)*TWO*C))
      RETURN
    2 A=-GC(IE)*SQRT((C+GC(IE)*CM)*(C-GC(IE)*CM+ONE)/
     *((C+ONE)*TWO*C))
      RETURN
      END
