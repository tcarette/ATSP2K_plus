*-----------------------------------------------------------------------
*        S A V E
*-----------------------------------------------------------------------
*
      SUBROUTINE SAVE(ICASE,C,K,I1,I2,I3,I4,JA,JB,IPTR)
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /FOUT/NOV(2),IOVLAP(10,2),NCOUNT(8),IFLAG,NIJ
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(8)

      IF (ICASE .LE. 2 .or. ICASE .EQ. 4 .or. ICASE .EQ. 5) THEN
*
*     Fk, Gk, L, or Z data
*
	 IF (I2 .GT. I4) THEN
	    II2 = I4
	    II4 = I2
	 ELSE
	    II2 = I2
	    II4 = I4
	 END IF
	 IPACK = (K*64 + II2)*64 + II4
	 IF (ICASE .NE. 4) THEN
	   WRITE(ISC(ICASE)) C,IPACK,JA,JB
	 ELSE
	   WRITE(ISC(ICASE)) C,IPACK,JA,JB,IPTR
	 END IF
      ELSE IF(ICASE .EQ. 3) THEN
	 J = 1
	 IMIN = I1
	 IF (I2 .LT. IMIN) THEN
	    IMIN=I2
	    J = 2
	 END IF
	 IF (I3 .LT. IMIN) THEN
	    IMIN = I3
	    J = 3
	 END IF
	 IF (I4 .LT. IMIN) THEN
	    IMIN = I4
	    J = 4
	 END IF
	 GO TO (10,20,30,40) J
10       II1 = I1
         II2 = I2
         II3 = I3
         II4 = I4
         Go to 50
	
20       II1 = I2
	 II2 = I1
	 II3 = I4
	 II4 = I3
	 GO TO 50

30       II1 = I3
	 II2 = I4
	 II3 = I1
	 II4 = I2
	 GO TO 50

40       II1 = I4
	 II2 = I3
	 II3 = I2
	 II4 = I1

50       IPACK = (((K*64+II1)*64+II2)*64+II3)*64+II4
	 WRITE(ISC(3)) C,IPACK,JA,JB,IPTR
      ELSE 
	 II1 = I1
	 II3 = I3
	 IF (I2 .GT. I4) THEN
	    II2 = I4
	    II4 = I2
	 ELSE
	    II2 = I2
	    II4 = I4
	 END IF
	 IF (ICASE .NE. 7) THEN
	    IF (I1 .GT. I3) THEN
	       II1 = I3
	       II3 = I1
	    END IF
	 END IF
*
* ... Because the k-value may be -1, for these integrals, a
*     value of k+1 is stored.
*
	 KK = K + 1
         IPACK = (((KK*64+II1)*64+II2)*64+II3)*64+II4
	 WRITE(ISC(ICASE)) C,IPACK,JA,JB
      END IF
      NCOUNT(ICASE) = NCOUNT(ICASE) + 1
      IFLAG = 1
      END
