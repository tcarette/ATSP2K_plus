*
*     -------------------------------------------------------------
*       C F G O 1
*     -------------------------------------------------------------
*
*     Modified by Gediminas Gaigalas,                September 1997
*
*
      SUBROUTINE CFGO1(NCFG,MAXORB,NCLOSD,QIAJCMP,QLJCOMP,QNJCOMP,
     :                 QNOC,QNELCSH,QNOCORB,QJ1,QIAJCLD,QLJCLSD)
*
*       Read configurations for one state, assuming orthogonality of
*       the orbitals
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=60,NWCD=20)
CAB      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
CAB     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
CAB      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
CAB     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
CAB     :       (QLJCLSD,LJCLSD(1))
      POINTER(QNOC,NOCCSH(*)),(QNELCSH,NELCSH(8,*)),
     :       (QNOCORB,NOCORB(8,*)),(QJ1,J1QNRD(15,*))
      POINTER(QIAJCMP,IAJCMP(*)),(QLJCOMP,LJCOMP(*)),
     :       (QNJCOMP,NJCOMP(*)),(QIAJCLD,IAJCLD(*)),
     :       (QLJCLSD,LJCLSD(*))
      CHARACTER*1 JAJCLD(3,NWCD),JAJCMP(3,NWD),JCQN(15),
     :            JQNG(15),J3QNG(15)
      CHARACTER EL(NWD)*3, LINE*72, HEAD*30
      DIMENSION J3QN(15),J2QN(15),J1QN(15)
*
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,JSC0,JSC(4)
*
    3 FORMAT(18(1X,A3))
    4 FORMAT(3A1)
**    5 FORMAT(8(1X,A3,1H(,I2,1H)))
    5 format(8(1X,A3,1X,I2,1X))
CGG    6 FORMAT(15(1X,I1,A1,I1))
   36 FORMAT(15(1X,A1,A1,A1))
    7 FORMAT(A72)
    8 FORMAT(A3)
    9 FORMAT(A15,F14.7)
   22 FORMAT(// 7H STATE ,' (WITH',I8,' CONFIGURATIONS):'/1H ,36(1H-)/)
   23 FORMAT(/10H THERE ARE,I3,21H ORBITALS AS FOLLOWS://
     : 5X,21(1X,A3):/5X,21(1X,A3))
   25 FORMAT(/14H CONFIGURATION,I5,' :'
     : ,8(1X,A3,1H(,I2,1H)))
CGG   26 FORMAT(4X,17H COUPLING SCHEME:,8(1X,4X,I1,A1,I1))
   26 FORMAT(4X,17H COUPLING SCHEME:,8(1X,4X,A1,A1,I1))
   46 FORMAT(4X,17H COUPLING SCHEME:,8(1X,4X,I1,A1,A1))
   27 FORMAT(32X,7(1X,4X,I1,A1,I1))
   47 FORMAT(32X,7(1X,4X,I1,A1,A1))
   28 FORMAT(/10H THERE ARE ,I3,31H CLOSED SUBSHELLS COMMON TO ALL ,
     :  27H CONFIGURATIONS AS FOLLOWS: //
     :  5X, 21(1X,A3))
*
* --- ANALYZE INPUT DATA
*
	  
	  CALL ANALY1(IREAD,IWRITE,NCLOSD,MAXORB,0,NCFG,EL)
*
* --- ALLOCATE MEMORY: NWFD = MAXORB
*
      NWFD = MAXORB
	  call alloc(qiajcmp,nwfd,4)
      call alloc(qljcomp,nwfd,4)
      call alloc(qnjcomp,nwfd,4)
      call alloc(qnoc,ncfg,4)
      call alloc(qnelcsh,8*ncfg,4)
      call alloc(qnocorb,8*ncfg,4)
      call alloc(qj1,15*ncfg,4)
      call alloc(qiajcld,nwcd,4)
      call alloc(qljclsd,nwcd,4)
*
      REWIND(IREAD)
*
* ---  Process the configuration data
*
      DO 30 I = 1,MAXORB
         READ(EL(I),8) IAJCMP(I)
         READ(EL(I),4) (JAJCMP(J,I),J=1,3)
30    CONTINUE
      WRITE(IWRITE,22) NCFG
      WRITE(IWRITE,23) MAXORB,(IAJCMP(I),I=1,MAXORB)
      DO 60 I=1,MAXORB
      IF (JAJCMP(1,I) .EQ. ' ') THEN
         JAJCMP(1,I) = JAJCMP(2,I)
         JAJCMP(2,I) = JAJCMP(3,I)
         JAJCMP(3,I) = ' '
      ENDIF
      LJCOMP(I) = LVAL(JAJCMP(2,I))
      NJCOMP(I) = ICHAR(JAJCMP(1,I)) - ICHAR('1') + 1
   60 CONTINUE
*
* --- READ HEADER CARD FOR THE CASE
*
      READ(IREAD,'(A30,I3,I4)') HEAD, NCI, NWI
*     WRITE(IOUT,7) HEAD
*
* --- READ IN THE COMMON SET OF CLOSED SUBSHELLS
*
      READ(IREAD,3) (EL(I),I=1,NCLOSD)
      DO 70 I=1,NCLOSD
         READ(EL(I),8) IAJCLD(I)
         READ(EL(I),4) (JAJCLD(J,I),J=1,3)
         J = 3
         IF (JAJCLD(1,I) .NE. ' ') J = 2
         LJCLSD(I) = LVAL(JAJCLD(J,I))
 70   CONTINUE
      WRITE(IWRITE,28) NCLOSD,(IAJCLD(I),I=1,NCLOSD)
*
* --- CHECK FORMAT OF CONFIGURATION FILE HEADER
*
 72   IF (NWI .GT. NCI) THEN
	READ(IREAD,'(A)')
	NWI = NWI - 20
	GO TO 72
      END IF
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
*     WRITE(IWRITE,25) I,(NOCORB(J,I),NELCSH(J,I),J=1,N)
      DO 61 J=1,N
      DO 61 JJ=1,MAXORB
   61 IF(NOCORB(J,I).EQ.IAJCMP(JJ)) NOCORB(J,I)=JJ
      M=2*N-1
      N1=N+1
CGG      READ(IREAD,6)    (J3QN(J),JCQN(J),J1QN(J),J=1,M)
      READ(IREAD,36)    (J3QNG(J),JCQN(J),JQNG(J),J=1,M)
*     WRITE(IWRITE,26) (J3QNG(J),JCQN(J),J1QN(J),J=1,N)
*     IF(N.GT.1) WRITE(IWRITE,27) (J3QN(J),JCQN(J),J1QN(J),J=N1,M)
      DO 62 J=1,M
      J2QN(J) = 2*LVAL(JCQN(J)) + 1
CGG
      J1QN(J) = NUMVAL(JQNG(J))
C      J3QN(J) = NUMVAL(J3QNG(J))
C      J3QN(J) = CHAR(J3QNG(J)+ICHAR('0'))
      J3QN(J) = ICHAR(J3QNG(J))-ICHAR('1') + 1
CGG
      J1QNRD(J,I) = (J3QN(J)*64 + J2QN(J))*64 + J1QN(J)
   62 CONTINUE
   63 CONTINUE
      CALL CFGTST(NCFG,QLJCOMP,QNOC,QNELCSH,QNOCORB,QJ1)
      RETURN
      END
