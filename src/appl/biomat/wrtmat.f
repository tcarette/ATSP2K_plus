*
*     ------------------------------------------------------------------
*	W R T M A T
*     ------------------------------------------------------------------
*
      SUBROUTINE WRTMAT(A,NROW,NCOL,NMROW,NMCOL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NMROW,NMCOL)
C
      DO 100 I=1,NROW
      WRITE(6,1010) I,(A(I,J),J=1,NCOL)
Cmrg 1010 FORMAT(1H0,I3,2X,4(1X,E14.8),/,(1H ,5X,4(1X,E14.8)))
 1010 FORMAT(1H0,I5,2X,4(1X,E14.8),/,(1H ,7X,4(1X,E14.8)))
  100 CONTINUE
      RETURN
      END

