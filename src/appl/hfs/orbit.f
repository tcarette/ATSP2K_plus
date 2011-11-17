*
*     -------------------------------------------------------------
*      O R B I T 
*     -------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      *
*     Vanderbilt University,  Nashville           September 1997   *
*
      SUBROUTINE ORBIT(L1,L2,I,A)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      A=ZERO
      IF((L1+L2).LT.1.OR.L1.NE.L2) RETURN
C      A=-DSQRT(DBLE(J1QN1(2*IHSH-1,2)*J1QN1(2*IHSH-1,3)))
C      A=A*DSQRT(DBLE(L1*(L1+1)*(L1+L1+1)))
      A=-DSQRT(DBLE(J1QN1(2*IHSH-1,2)))
      A=A*DSQRT(TWO*DBLE(L1*(L1+1)*(L1+L1+1)))
      RETURN
      END
