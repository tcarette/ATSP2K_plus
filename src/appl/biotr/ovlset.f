*
*-----------------------------------------------------------------------
*        O V L S E T
*-----------------------------------------------------------------------
*
      SUBROUTINE OVLSET(N,I1,I2,IOV,I3,I4,JOV,IPTR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,
     : IALL,JSC(3),ISCW
 
      write(Iscw,*) ' Overlap Integrals have been encountered'
      write(Iscw,*) ' which are not allowed in this implementation'
      STOP
      END
