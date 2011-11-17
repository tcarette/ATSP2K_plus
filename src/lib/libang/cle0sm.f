*
*     -------------------------------------------------------------
*      C L E 0 S M
*     -------------------------------------------------------------
*
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
*                                                 ---         ---  *
*                                                 I  Q   S  C   I  *
*     CLEBSCH - GORDAN COEFFICIENT:               I             I  *
*                                                 I  QM  0  CM  I  *
*                                                 ---         ---  *
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vilnius,  Lithuania                          December 1993   *
*
      SUBROUTINE CLE0SM(Q,QM,S,C,CM,A)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      A=ZERO
      IIQ=TWO*Q+TENTH
      IIC=TWO*C+TENTH
      IIS=TWO*S+TENTH
      IF(ITTK(IIQ,IIC,IIS).EQ.0)RETURN
      IF(S.LT.EPS)GO TO 1
      CALL C1E0SM(Q,QM,C,CM,A)
      RETURN
    1 IF((Q+TENTH).LT.ABS(QM))RETURN
      IF((C+TENTH).LT.ABS(CM))RETURN
      IF(ABS(Q-C).GT.EPS)RETURN
      IF(ABS(QM-CM).GT.EPS)RETURN
      A=ONE
      RETURN
      END
