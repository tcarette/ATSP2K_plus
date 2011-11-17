*
*     ------------------------------------------------------------------
*	I E L S U M
*     ------------------------------------------------------------------
*
      FUNCTION IELSUM(IVEC,NELMNT)
*
* Sum elements of integer array
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IVEC(*)
*
      ISUM = 0
C?    WRITE(6,*) ' IELSUM, NELMNT ', NELMNT
      DO 100 IEL = 1, NELMNT
       ISUM = ISUM + IVEC(IEL)
C?     write(6,*) ' IELSUM IEL IVEC ISUM '
C?     write(6,*) IEL,IVEC(IEL),ISUM
  100 CONTINUE
*
      IELSUM = ISUM
*
      RETURN
      END
