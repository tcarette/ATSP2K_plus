*
*     ------------------------------------------------------------------
*	G T R A C 1
*     ------------------------------------------------------------------
*
      SUBROUTINE GTRAC1(I,L,IFIRST,NFOUND,BUF,
     &           IBUF,LBUF,JLEI,LU,NTESTG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
* Obtain racah coefficient for given I and l starting 
* from element IFIRST
*
* JLEI = 0 => include all terms with given I
* JLEI = 1 => include only terms with J<I 
*
*
C     write(6,*) ' Old GTRACX '
C     CALL  GTRACX(I,L,IFIRST,NFOUND,BUF,
C    &                IBUF,LBUF,JLEI,LU,NTESTG)
*
C     write(6,*) ' new gtracx '
C     CALL  GTRACXN(I,L,IFIRST,NFOUND,BUF,
C    &                IBUF,LBUF,JLEI,LU,NTESTG)
* For new version where the coupling coefficients reside in core,
* the program must know which of the two lists that should be read.
* This is done by assuming that LU is as ' in the bad old days'
* i.e
*         LU = 15 => Left (initial state)
*         LU = 16 => right (final state)
      IF(LU.EQ.15) THEN
        ILIST = 1
      ELSE IF (LU.EQ.16) THEN
        ILIST = 2
      ELSE
        WRITE(6,*) ' PROBLEM in GTRAC1, Unrecognized LU = ',LU
        STOP'GTRAC1: Undefined LU'
      END IF
*
      CALL  GTRACXVN(I,L,IFIRST,NFOUND,BUF,
     &                IBUF,LBUF,JLEI,ILIST,NTESTG)
*
      RETURN
      END 
