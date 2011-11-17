*
*
*     ------------------------------------------------------------------
*	S C A L V E
*     ------------------------------------------------------------------
*
      SUBROUTINE SCALVE(VECTOR,FACTOR,NDIM)
C
C CALCULATE SCALAR(FACTOR) TIMES VECTOR
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VECTOR(1)
C
      DO 100 I=1,NDIM
       VECTOR(I)=VECTOR(I)*FACTOR
  100 CONTINUE
C
      RETURN
      END
