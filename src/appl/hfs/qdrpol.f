*
*     -------------------------------------------------------------
*      Q D R P O L 
*     -------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      *
*     Vanderbilt University,  Nashville           September 1997   *
*
      SUBROUTINE QDRPOL(L1,L2,I,A)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      A=ZERO
      IF(ITTK(L1,L2,2).EQ.0) RETURN
C      A=-DSQRT(DBLE(J1QN1(2*IHSH-1,2)*J1QN1(2*IHSH-1,3)))
      A=-DSQRT(TWO*DBLE(J1QN1(2*IHSH-1,2)))
      A=A*RME(L1,L2,2)
      IF(MOD(L2-L1+2,4).NE.0) A=-A
      RETURN
      END
