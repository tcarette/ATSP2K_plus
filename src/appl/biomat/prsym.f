*
*     ------------------------------------------------------------------
*	P R S Y M   
*     ------------------------------------------------------------------
*
      SUBROUTINE PRSYM(A,MATDIM)
C PRINT LOWER HALF OF A SYMMETRIC MATRIX OF DIMENSION MATDIM.
C THE LOWER HALF OF THE MATRIX IS SUPPOSED TO BE IN VECTOR A.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1)
      JSTART=1
      JSTOP=0
      DO 100 I=1,MATDIM
        JSTART=JSTART+I-1
        JSTOP=JSTOP +I
        WRITE(6,1010) I,(A(J),J=JSTART,JSTOP)
  100 CONTINUE
      RETURN
 1010 FORMAT(1H0,2X,I3,5(1X,E13.7),/,(1H ,5X,5(1X,E13.7)))
      END
