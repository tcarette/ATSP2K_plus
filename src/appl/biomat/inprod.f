*
*     ------------------------------------------------------------------
*	I N P R O D
*     ------------------------------------------------------------------
*
      DOUBLE PRECISION FUNCTION INPROD(A,B,NDIM)
C      CALCULATE SCALAR PRODUCT BETWEEN TO VECTORS A,B
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(2),B(2)
C
      INPROD=0.0D0
      DO 100 I=1,NDIM
       INPROD=INPROD+A(I)*B(I)
  100 CONTINUE
C
      RETURN
      END
