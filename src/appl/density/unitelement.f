*
*     -------------------------------------------------------------
*      U N I T E L E M E N T
*     -------------------------------------------------------------
*
*
      SUBROUTINE UNITELEMENT(L1,L2,I,A)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)

      A=-DSQRT(DBLE(4*L1+2))
      RETURN
      END
