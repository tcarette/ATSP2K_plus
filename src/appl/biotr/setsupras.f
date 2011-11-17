*
*     --------------------------------------------------------------
*	S E T S U P R A S 
*     --------------------------------------------------------------
*
* --- Sets up the RAS information of shells. 
*
*     Written by G. Gaigalas,                                      *
*     Vanderbilt University,  Nashville           September 1997   *
*
      SUBROUTINE SETSUPRAS(IK,NCLOSD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=128)
      CHARACTER*3 elras,elrasi,elrasf
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc(2),iuw(2),iul(2),iuj(2),iut(2)
      COMMON /RAS/ ninac(11,2),nras1(11,10,2),nras2(11,2),
     :             nras3(11,10,2),NGAS1(2),NGAS3(2),
     :             nl(11,2),elras(2*nwd),elrasi(nwd),elrasf(nwd),
     :             itab(nwd,2),lmax(2)
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      cOMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      LM=0
      DO 1 I=1,11
        NINAC(J,IK)=0
        NRAS2(J,IK)=0
    1 CONTINUE
      IF(NCLOSD.NE.0) THEN
        DO 2 I=1,NCLOSD
          J=LJCLSD(I)+1
          IF(LM.LT.J) LM=J
	  NINAC(J,IK)=NINAC(J,IK)+1
    2   CONTINUE
      ENDIF
      IF(MAXORB.NE.0) THEN
        DO 3 I=1,MAXORB
          J=LJCOMP(I)+1
          IF(LM.LT.J) LM=J
	  NRAS2(J,IK)=NRAS2(J,IK)+1
    3   CONTINUE
      ENDIF
      IF(LM.GT.11) THEN
	WRITE(ISCW,*) ' l-value too large in RASIN (SETSUPRAS)  '
	STOP
      ELSE
        LMAX(IK)=LM
      ENDIF
      RETURN 
      END
