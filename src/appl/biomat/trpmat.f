*
*     ------------------------------------------------------------------
*	T R P M A T
*     ------------------------------------------------------------------
*
      SUBROUTINE TRPMAT(XIN,NROW,NCOL,XOUT)
C
C XOUT(I,J) = XIN(J,I)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XIN(NROW,NCOL),XOUT(NCOL,NROW)
C
      DO 200 IROW =1, NROW
        DO 100 ICOL = 1, NCOL
          XOUT(ICOL,IROW) = XIN(IROW,ICOL)
  100   CONTINUE
  200 CONTINUE
C
      RETURN
      END
