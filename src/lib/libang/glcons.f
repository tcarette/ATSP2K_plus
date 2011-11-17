*
*     ------------------------------------------------------------------
*       B L O C K   D A T A    G L C O N S
*     ------------------------------------------------------------------
*
      BLOCK DATA GLCONS
C
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
C
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
C
C     SET GLOBAL REAL CONSTANTS
C
      DATA ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS/
     :     0.0D 00,0.1D 00,0.5D 00,
     :     1.0D 00,2.0D 00,3.0D 00,
     :     4.0D 00,
     :     7.0D 00,1.1D 01,1.0D-08/
C
      END

