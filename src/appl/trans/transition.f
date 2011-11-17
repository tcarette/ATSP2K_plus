*
*     -------------------------------------------------------------
*      T R A N S I T I O N 
*     -------------------------------------------------------------
*
*     CALCULATE TRANSITION OPERATORS                               *
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vanderbilt University,  Nashville           September 1997   *
*
      SUBROUTINE TRANSITION(L1,L2,I,A)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      LOGICAL REL,VOK
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON /EMS/IEM(4),IFL,JI,JF,LAM,REL,VOK
      IF(IFL.NE.4) THEN
        A=-DSQRT(DBLE(TWO*J1QN1(2*IHSH-1,2)))
      ELSE
        A=-DSQRT(DBLE(J1QN1(2*IHSH-1,2)*J1QN1(2*IHSH-1,3)))
      ENDIF
      RETURN
      END
