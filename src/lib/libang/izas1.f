*
*     -------------------------------------------------------------
*      I Z A S 1
*     -------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      *
*     Vilnius,  Lithuania                             March 1995   *
*
      FUNCTION IZAS1(IB,QB,IK,QK)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      IZAS1=0
      IQB=TWO*ABS(QB)+TENTH
      IF(IQB.GT.IB)RETURN
      IF(MOD(IB+IQB,2).NE.0)RETURN
      IQK=TWO*ABS(QK)+TENTH
      IF(IQK.GT.IK)RETURN
      IF(MOD(IK+IQK,2).NE.0)RETURN
      IZAS1=1
      RETURN
      END
