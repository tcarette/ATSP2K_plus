*
*     -------------------------------------------------------------
*      C 0 T 5 S
*     -------------------------------------------------------------
*
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
*                                                ---          ---  *
*                                                I  Q  1/2  C   I  *
*     CLEBSCH - GORDAN COEFFICIENT:              I              I  *
*                                                I  QM  SM  CM  I  *
*                                                ---          ---  *
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vilnius,  Lithuania                          December 1993   *
*
      SUBROUTINE C0T5S(Q,QM,SM,C,CM,A)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DIMENSION GC(2)
      GC(1)=ONE
      GC(2)=-ONE
      A=ZERO
      IIQ=TWO*Q+TENTH
      IIC=TWO*C+TENTH
      IF(ITTK(IIQ,IIC,1).EQ.0)RETURN
      IF(DABS(QM+SM-CM).GT.EPS)RETURN
      IF((HALF+TENTH).LT.DABS(SM))RETURN
      IF((Q+TENTH).LT.DABS(QM))RETURN
      IF((C+TENTH).LT.DABS(CM))RETURN
      IE=DABS(HALF-SM)+ONE+TENTH
      IF(DABS(Q+HALF-C).LT.EPS) THEN
        A=DSQRT((C+GC(IE)*CM)/(TWO*C))
      ELSE
        IF(DABS(Q-HALF-C).GT.EPS)RETURN
        A=-GC(IE)*DSQRT((C-GC(IE)*CM+ONE)/(TWO*C+TWO))
      ENDIF
      RETURN
      END
