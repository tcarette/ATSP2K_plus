*
*
*     ------------------------------------------------------------------
*	I N V M A T
*     ------------------------------------------------------------------
*
      SUBROUTINE INVMAT(A,B,MATDIM,NDIM)
C FIND INVERSE OF MATRIX A
C INPUT :
C        A : MATRIX TO BE INVERTED
C        B : SCRATCH ARRAY
C        MATDIM : PHYSICAL DIMENSION OF MATRICES
C        NDIM :   DIMENSION OF SUBMATRIX TO BE INVERTED
C
C OUTPUT : A : INVERSE MATRIX ( ORIGINAL MATRIX THUS DESTROYED )
C WARNINGS ARE ISSUED IN CASE OF CONVERGENCE PROBLEMS )
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(MATDIM,MATDIM),B(MATDIM,MATDIM)
C
      ITEST=0
      IF(NDIM.EQ.1)THEN
        IF(A(1,1) .NE. 0.0D0 ) THEN
           A(1,1) = 1.0D0/A(1,1)
        ELSE
           ITEST = 1
        END IF
      ELSE
        DETERM=0.0D0
        EPSIL=0.0D0
        CALL BNDINV(A,B,NDIM,DETERM,EPSIL,ITEST,MATDIM)
      END IF
C
      IF( ITEST .NE. 0 ) THEN
        WRITE (6,'(A,I3)') ' INVERSION PROBLEM NUMBER..',ITEST
      END IF
      NTEST = 1   
      IF ( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' INVERTED MATRIX '
        CALL WRTMAT(A,NDIM,NDIM,MATDIM,MATDIM)
      END IF
C
      RETURN
      END
