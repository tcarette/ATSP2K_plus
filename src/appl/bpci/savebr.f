*------------------------------------------------------------------------
*        S A V E B R
*------------------------------------------------------------------------
*
*     Three matrices are created -- LS (non-relativistic or with
*     relativistic shift), LSJ1 (Breit-Pauli operators other
*     than spin-spin), LSJ2 (spin-spin)
*
*     An LS calculation only uses the first whereas LSJ creates
*     J-dependent matrices from all three.
*
      SUBROUTINE SAVE(ICASE,C,K,I1,I2,I3,I4,JA,JB,IPTR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      POINTER (qh,h(ncfg,3)),(qjan,jan(ncfg,3))
      COMMON /buffer/qh, qjan, nrow(3),iflag(3)
      COMMON /INFORM/IREAD,IWRITE,IOUT,ISC(8),ISCW
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      DOUBLE PRECISION intval
*
*     print *,'Save: icase,c,k,i1,i2,i3,i4,ja,jb,iptr',
*    :        icase,c,k,i1,i2,i3,i4,ja,jb,iptr
      if (icase .le. 4) then
	lcase = 1
      else if (icase .le. 7) then
	lcase = 2
      ELSEIF(ICASE.EQ.9) THEN
	LCASE=1
      else
	lcase = 3
      end if
*
*     .. obtain canonical form
      IF (icase .LE. 2 .or. icase .EQ. 4 .or. icase .eq. 5) THEN
	 IF (I2 .GT. I4) THEN
	    II2 = I4
	    II4 = I2
	 ELSE
	    II2 = I2
	    II4 = I4
	 END IF
	 
      ELSE if (icase .eq. 3) then
*
*     Rk data
*
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
      END IF
50    v= intval(ICASE,K,II1,II2,II3,II4)
      if (iflag(lcase) .eq.0) then
*       .. this is the first for this matrix element
	nrow(lcase) = nrow(lcase) + 1
	jan(nrow(lcase),lcase) = ja
        iflag(lcase) = 1
      end if
      h(nrow(lcase),lcase) = h(nrow(lcase),lcase) + c*v
      END
