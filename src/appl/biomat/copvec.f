*
*     ------------------------------------------------------------------
*	C O P V E C
*     ------------------------------------------------------------------
*
      SUBROUTINE COPVEC(FROM,TO,NDIM)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION FROM(ndim),TO(ndim)
C
      DO 100 I=1,NDIM
       TO(I)=FROM(I)
  100 CONTINUE
C
      RETURN
      END
