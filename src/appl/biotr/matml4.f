*
*     ------------------------------------------------------------------
*	M A T M L 4
*     ------------------------------------------------------------------
*
      SUBROUTINE MATML4(C,A,B,NCROW,NCCOL,NAROW,NACOL,
     &                  NBROW,NBCOL,ITRNSP )
C
C MULTIPLY A AND B TO GIVE C
C
C     C = A * B             FOR ITRNSP = 0
C
C     C = A(TRANSPOSED) * B FOR ITRNSP = 1
C
C     C = A * B(TRANSPOSED) FOR ITRNSP = 2
C
C... JEPPE OLSEN, LAST REVISION JULY 24 1987
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NAROW,NACOL),B(NBROW,NBCOL)
      DIMENSION C(NCROW,NCCOL)
C
      NTEST = 0
      IF ( NTEST .NE. 0 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' A AND B MATRIX FROM MATML4 '
        WRITE(6,*)
        CALL WRTMAT(A,NAROW,NACOL,NAROW,NACOL)
        CALL WRTMAT(B,NBROW,NBCOL,NBROW,NBCOL)
        WRITE(6,*)      ' NCROW NCCOL NAROW NACOL NBROW NBCOL '
        WRITE(6,'(6I6)')  NCROW,NCCOL,NAROW,NACOL,NBROW,NBCOL
      END IF
C
      CALL SETVEC(C,0.0D0,NCROW*NCCOL)
C
      IF( ITRNSP .NE. 0 ) GOTO 001
        DO 50 J = 1,NCCOL
          DO 40 K = 1,NBROW
            BKJ = B(K,J)
            DO 30 I = 1, NCROW
              C(I,J) = C(I,J) + A(I,K)*BKJ
  30        CONTINUE
  40      CONTINUE
  50    CONTINUE
C
C
  001 CONTINUE
C
      IF ( ITRNSP .NE. 1 ) GOTO 101
C... C = A(T) * B
         DO 150 J = 1, NCCOL
           DO 140 K = 1, NBROW
             BKJ = B(K,J)
             DO 130 I = 1, NCROW
               C(I,J) = C(I,J) + A(K,I)*BKJ
  130        CONTINUE
  140      CONTINUE
  150    CONTINUE
C
  101 CONTINUE
C
      IF ( ITRNSP .NE. 2 ) GOTO 201
C... C = A*B(T)
        DO 250 J = 1,NCCOL
          DO 240 K = 1,NBCOL
            BJK = B(J,K)
            DO 230 I = 1, NCROW
              C(I,J) = C(I,J) + A(I,K)*BJK
 230        CONTINUE
 240      CONTINUE
 250    CONTINUE
C
C
  201 CONTINUE
C
      IF ( NTEST .NE. 0 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' C MATRIX FROM MATML4 '
        WRITE(6,*)
        CALL WRTMAT(C,NCROW,NCCOL,NCROW,NCCOL)
      END IF
C
      RETURN
      END
