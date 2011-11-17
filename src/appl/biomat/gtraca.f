*
*     ------------------------------------------------------------------
*	G T R A C A
*     ------------------------------------------------------------------
*
      SUBROUTINE GTRACA(I,L,IFIRST,NFOUND,BUF,
     &           IBUF,LBUF,JLEI,LU)
*
* 1s inactive, 2 electrons in (1s 1p), Singlet S
*
* The following is the list in the order : l  j i   r l coef
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NCOUP =2)
      DIMENSION IRACAH(5,NCOUP)
      DIMENSION  RACAH(NCOUP)
      DATA  IRACAH/
     &0, 2,2, 1,1,
     &1, 1,1, 2,2/
*234567
*. Output
      DIMENSION BUF(LBUF),IBUF(4,LBUF)
      PARAMETER (SQRT2 = 1.414D0     )
*
      DATA RACAH/
     & 2.0D0,2.0D0/
*
*. Position so IFIRST -1 elements have been skipped
*
     
      IF(IFIRST.NE.1) THEN
       ISKIP = 0
       DO 100 IELMNT = 1, NCOUP
        JELMNT = IELMNT
        IF(IRACAH(3,IELMNT) .EQ. I .AND.
     &     IRACAH(1,IELMNT).EQ.L ) THEN
*
          IF(JLEI.EQ.0) THEN
            ISKIP = ISKIP +1
          ELSE IF (JLEI.EQ.1.AND.
     &    IRACAH(2,IELMNT).LT.IRACAH(3,IELMNT))THEN
            ISKIP = ISKIP + 1
          ELSE IF (JLEI.EQ.2.AND.
     &    IRACAH(2,IELMNT).EQ.IRACAH(3,IELMNT))THEN
            ISKIP = ISKIP + 1
          ELSE IF (JLEI.EQ.3.AND.
     &    IRACAH(2,IELMNT).LE.IRACAH(3,IELMNT))THEN
            ISKIP = ISKIP + 1
          END IF
          IF(ISKIP.EQ.IFIRST-1) GOTO 101
*
        END IF
  100  CONTINUE
  101  CONTINUE
      ELSE IF ( IFIRST .EQ. 1 ) THEN
         JELMNT = 0
      END IF
*
* Obtain the next elements, atmost LBUF
*
       NFOUND = 0
       DO 200 IELMNT = JELMNT+1, NCOUP
        IF(IRACAH(3,IELMNT) .EQ. I .AND.
     &     IRACAH(1,IELMNT).EQ.L ) THEN
          IF(JLEI.EQ.0) THEN
           NFOUND = NFOUND + 1
           IBUF(1,NFOUND) = IRACAH(2,IELMNT)
           IBUF(2,NFOUND) = IRACAH(3,IELMNT)
           IBUF(3,NFOUND) = IRACAH(4,IELMNT)
           IBUF(4,NFOUND) = IRACAH(5,IELMNT)
            BUF(  NFOUND) =  RACAH(  IELMNT)
          ELSE IF (JLEI.EQ.1.AND.
     &    IRACAH(2,IELMNT).LT.IRACAH(3,IELMNT))THEN
           NFOUND = NFOUND + 1
           IBUF(1,NFOUND) = IRACAH(2,IELMNT)
           IBUF(2,NFOUND) = IRACAH(3,IELMNT)
           IBUF(3,NFOUND) = IRACAH(4,IELMNT)
           IBUF(4,NFOUND) = IRACAH(5,IELMNT)
            BUF(  NFOUND) =  RACAH(  IELMNT)
          ELSE IF (JLEI.EQ.2.AND.
     &    IRACAH(2,IELMNT).EQ.IRACAH(3,IELMNT))THEN
           NFOUND = NFOUND + 1
           IBUF(1,NFOUND) = IRACAH(2,IELMNT)
           IBUF(2,NFOUND) = IRACAH(3,IELMNT)
           IBUF(3,NFOUND) = IRACAH(4,IELMNT)
           IBUF(4,NFOUND) = IRACAH(5,IELMNT)
            BUF(  NFOUND) =  RACAH(  IELMNT)
          ELSE IF (JLEI.EQ.3.AND.
     &    IRACAH(2,IELMNT).LE.IRACAH(3,IELMNT))THEN
           NFOUND = NFOUND + 1
           IBUF(1,NFOUND) = IRACAH(2,IELMNT)
           IBUF(2,NFOUND) = IRACAH(3,IELMNT)
           IBUF(3,NFOUND) = IRACAH(4,IELMNT)
           IBUF(4,NFOUND) = IRACAH(5,IELMNT)
            BUF(  NFOUND) =  RACAH(  IELMNT)
          END IF
          IF(NFOUND.EQ.LBUF   ) GOTO 201
        END IF
  200  CONTINUE
  201  CONTINUE
*
       NTEST = 0
       IF( NTEST.NE.0) THEN
        WRITE(6,*) ' GTRAC1 to your service '
        WRITE(6,*) ' ======================='
        WRITE(6,*) 'Number of elements obtained',NFOUND
        WRITE(6,*) ' IBUF and BUF '
        DO 300 IELMNT = 1, NFOUND
          WRITE(6,'(E12.7,3X,4I4)') 
     &    BUF(IELMNT),(IBUF(K,IELMNT),K=1,4)
  300  CONTINUE
*
      END IF
      RETURN
      END
