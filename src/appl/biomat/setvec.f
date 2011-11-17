*
*     ------------------------------------------------------------------
*	S E T V E C
*     ------------------------------------------------------------------
*
      SUBROUTINE SETVEC(VECTOR,VALUE,NDIM)
C
C VECTOR (*) = VALUE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VECTOR(*)
C
      DO 10 I=1,NDIM
   10 VECTOR(I) = VALUE
C
      RETURN
      END
