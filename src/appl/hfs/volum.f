*
*     -------------------------------------------------------------
*      V O L U M 
*     -------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      *
*     Vanderbilt University,  Nashville           September 1997   *
*
      SUBROUTINE VOLUM(L1,L2,I,A)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      A=ZERO
      IF(L1.NE.0.OR.L2.NE.0) RETURN
      A=-DSQRT(DBLE(J1QN1(2*IHSH-1,2)*J1QN1(2*IHSH-1,3)))
      RETURN
      END
