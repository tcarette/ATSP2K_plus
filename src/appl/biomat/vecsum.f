
*
*     ------------------------------------------------------------------
*	V E C S U M
*     ------------------------------------------------------------------
*
      SUBROUTINE VECSUM(C,A,B,FACA,FACB,NDIM)
C
C     CACLULATE THE VECTOR C(I)=FACA*A(I)+FACB*B(I)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1   ),B(1   ),C(1   )
*
      IF(FACA.NE.0.0D0.AND.FACB.NE.0.0D0) THEN
        DO 100 I=1,NDIM
          S=FACA*A(I)+FACB*B(I)
          C(I)=S
  100   CONTINUE
*
      ELSE IF(FACA.EQ.0.0D0.AND.FACB.NE.0.0D0) THEN
        DO 200 I=1,NDIM
          S=FACB*B(I)
          C(I)=S
  200   CONTINUE
*
      ELSE IF(FACA.NE.0.0D0.AND.FACB.EQ.0.0D0) THEN
        DO 300 I=1,NDIM
          S=FACA*A(I)
          C(I)=S
  300   CONTINUE
*
      ELSE IF(FACA.EQ.0.0D0.AND.FACB.EQ.0.0D0) THEN
        DO 400 I=1,NDIM
          C(I)=0.0D0
  400   CONTINUE
      END IF
C
      RETURN
      END
