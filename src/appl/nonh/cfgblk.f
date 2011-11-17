*     ------------------------------------------------------------------
*       C F G B L K
*     ------------------------------------------------------------------
*
      SUBROUTINE CFGBLK(ncfg,maxorb,QIAJCMP,QLJCOMP,QNJCOMP,QNOC,
     :                  QNELCSH,QNOCORB,QJ1,QIAJCLD,QLJCLSD,term)
*
*       Read configurations for one block 
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*     IMPLICIT INTEGER (Q)
      PARAMETER (NWD=70,NWCD=20)
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      CHARACTER EL(NWD)*3, LINE*72, HEAD*30
      DIMENSION J3QN(15),J2QN(15),J1QN(15)
      CHARACTER*1 JAJCLD(3,NWCD),JAJCMP(3,NWD),JCQN(15),
     :            JQNG(15),J3QNG(15)
*
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,JSC0,JSC(4)
      COMMON /CLOSED/B1ELC(4),NCLOSD,IBK
     : IALL,ISCW, state
      character term*3
      character*1 RD1, RD2, parity
*
    5 FORMAT(8(1X,A3,1X,I2,1X))
    6 FORMAT(15(1X,A1,A1,A1))
    7 FORMAT(A72)
*
*     .. skip the first two lines
      read(iread,7) 
      read(iread,7) 
*
* --- READ IN (AND PRINT OUT) CONFIGURATIONS ETC. FOR THE STATE UNDER
* --- CONSIDERATION
*
      DO 63 I=1,NCFG
      READ(IREAD,7) LINE
      N=0
      J=2
   65 IF (LINE(J:J+2) .NE. '   ' .and. N .LT. (8)) THEN
         N = N + 1
         J = J +8
         GO TO 65
      END IF
      NOCCSH(I) = N
      READ(LINE,5)       (NOCORB(J,I),NELCSH(J,I),J=1,N)
      DO 61 J=1,N
      DO 61 JJ=1,MAXORB
   61 IF(NOCORB(J,I).EQ.IAJCMP(JJ)) NOCORB(J,I)=JJ
      M=2*N-1
      N1=N+1
      read(iread,7) line
      READ(line,6)    (J3QNG(J),JCQN(J),JQNG(J),J=1,M)

      If (i .eq. 1) then
*       Determine parity
        ip = 0
        Do iip = 1,n
          iorb = nocorb(iip,1)
          iq   = nelcsh(iip,1)
          ip = ip + ljcomp(iorb)*iq
*         print *, iorb, iq, ip, ljcomp(iorb)
        end do
        if ((ip/2)*2 .eq. ip) then
          parity ='e'
        else
          parity ='o'
        end if
*       In small cases, term might be three characters as in 1D2.
*       In such cases, parity replaces the third character.
        len = len_trim(line)
        term = adjustl(line(len-2:len))
        term(3:3) = parity
	 print *, 'processing ', term, ' with ',ncfg, 'configurations'
      end if


      DO 62 J=1,M
      J1QN(J) = NUMVAL(JQNG(J))
      J2QN(J) = 2*LVAL(JCQN(J)) + 1
      J3QN(J) = ICHAR(J3QNG(J))-ICHAR('1') + 1
      J1QNRD(J,I) = (J3QN(J)*64 + J2QN(J))*64 + J1QN(J)
   62 CONTINUE
   63 CONTINUE
*     .. skip * line
      read(iread,7)
*
*     . check the coupling
      CALL CFGTST(NCFG,QLJCOMP,QNOC,QNELCSH,QNOCORB,QJ1)
*
      RETURN
      END

