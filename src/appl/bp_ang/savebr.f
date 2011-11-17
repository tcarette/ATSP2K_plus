*------------------------------------------------------------------------
*        S A V E B R
*------------------------------------------------------------------------
*
      SUBROUTINE SAVE(ICASE,C,K,I1,I2,I3,I4,JA,JB,IPTR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER (LSDIM=30000)
      POINTER (qcn,cn(lsdim)),(qinptr,inptr(lsdim)),
     :        (qnijptr,nijptr(lsdim)),(qjan,jan(lsdim)),
     :        (qjbn,jbn(lsdim)),(qintptr,intptr(0:2*lmax+1,7)),
     :        (qpackn,ipackn(1)),(qlused,lused(1)),(qico,ico(1))
      COMMON /buffer/qcn,qinptr,qpackn,qlused,qintptr,lmax,qnijptr,
     :               qjan,qjbn,qico
      POINTER  (qjptr, jptr(1))
      COMMON /fout/n,ntot,iflag,nih,nij,qjptr
      !COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(4),IALL,JSC(3),ISCW
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(8),ISCW

      LOGICAL lused
*
!      write(99,*)'SAVE:', ICASE,C,K,I1,I2,I3,I4,JA,JB,IPTR
      if (n .eq. LSDIM) then
*        .. write data to disk
         new = n
         write(50) new,(cn(j),j=1,new),(inptr(j),j=1,new)
         !write(*,*) new,(cn(j),j=1,new),(inptr(j),j=1,new)
         n = 1
      else
         n = n + 1
      end if

*     nij is the number of the current coefficient from the beginning
*     n is the number withing the current block

      nij = nij +1

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
50    CONTINUE
      klocal = k
*      .. pack the data
       if (icase.le.2) then
         int = (klocal*64 + ii2)*64 + ii4
       else if (icase.eq.4 .or. icase.eq.5) then
         int = ii2*64 + ii4
       else
         if (icase.eq.6.or.icase.eq.8.or.icase.eq.9) klocal = klocal+1
         int = (((klocal*64+II1)*64+II2)*64+II3)*64+II4
       endif
       IFLAG = 1
       Cn(n) = C
!       write(99,*) 'Isearch:', icase,int, klocal,ii1,ii2,ii3,ii4
       INPTR(n) = isearch(icase,int,qpackn,qintptr,lmax)
       if (icase .eq. 8) inptr(n) = -inptr(n)
      END
