*
*     ------------------------------------------------------------------
*       O R T H O G G
*     ------------------------------------------------------------------
*
      SUBROUTINE ORTHOGG(LET)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(4),IALL,JSC(3),ISCW
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/SIGNF /SIGNFA
*
*     THIS SUBROUTINE CHECKS FOR POSSIBLE ORTHOGONALITY DUE TO
*     COUPLING DIFFERENCES OR UNEVEN PARITY
*
  102 FORMAT(52H ORTHOGONALITY IN COUPLING SCHEMES OF CONFIGURATIONS)
  103 FORMAT(59H THE TWO CONFIGURATIONS HAVE DIFFERING NUMBERS OF ELECTR
     :ONS)
  104 FORMAT(51H THE TWO CONFIGURATIONS HAVE DIFFERING TOTAL PARITY)
*
* --- DO PSI AND PSIP CONTAIN THE SAME NUMBERS OF ELECTRONS
*     DO PSI AND PSIP HAVE THE SAME TOTAL PARITY
*
      N5=0
      N6=0
      N7=0
      DO 1 I=1,IHSH
        L1=LJ(I)
        L2=NOSH1(I)
        L3=NOSH2(I)
        N5=N5+L2
        N6=N6+L3
        N7=N7+L1*(L2-L3)
    1 CONTINUE
*
*     CHECK ON NUMBER OF ELECTRONS
*
      IF ((N5-N6).NE.0) THEN
        WRITE(IWRITE,103)
        LET=0
*
*     CHECK ON PARITY
*
      ELSE
        IF((N7-N7/2*2).NE.0) THEN
          WRITE(IWRITE,104)
          LET=0
        ELSE
          N72=N7/2
          SIGNFA=1.D0
          IF( (N72-(N72/2)*2).NE.0 ) SIGNFA=-SIGNFA
*
* --- COUPLING ORTHOGONALITY TEST FOR FIRST TWO SHELLS
*
          LET=1
        ENDIF
      ENDIF
      RETURN
      END
*
*     ------------------------------------------------------------------
*       O R T H O G G B
*     ------------------------------------------------------------------
*
      SUBROUTINE ORTHOGGB(LET,INCL)
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(8),ISCW
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON/IMAGNT/CONST,CONSOO,CONSS,ISPORB,ISOORB,ISPSPN,
     :     IREL,ISTRICT,IZOUT,IELST,ITENPR
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      LOGICAL INCL
*
*     THIS SUBROUTINE CHECKS FOR POSSIBLE ORTHOGONALITY DUE TO
*     COUPLING DIFFERENCES OR UNEVEN PARITY
*
  101 FORMAT(37H DIFFERING RESULTANT ANGULAR MOMENTUM)
  102 FORMAT(52H ORTHOGONALITY IN COUPLING SCHEMES OF CONFIGURATIONS)
  103 FORMAT(59H THE TWO CONFIGURATIONS HAVE DIFFERING NUMBERS OF ELECTR
     :ONS)
  104 FORMAT(51H THE TWO CONFIGURATIONS HAVE DIFFERING TOTAL PARITY)
*
* --- DO PSI AND PSIP CONTAIN THE SAME NUMBERS OF ELECTRONS
*     DO PSI AND PSIP HAVE THE SAME TOTAL PARITY
*
      N5=0
      N6=0
      N7=0
      IELST=1
      DO 20 I=1,IHSH
      L1=LJ(I)
      L2=NOSH1(I)
      L3=NOSH2(I)
      N5=N5+L2
      N6=N6+L3
      N7=N7+L1*(L2-L3)
   20 CONTINUE
*
*     CHECK ON NUMBER OF ELECTRONS
*
      IF (N5-N6) 21,22,21
   21 IF(IBUG2-1) 11,28,28
   28 WRITE(IWRITE,103)
      GO TO 11
*
*     CHECK ON PARITY
*
   22 IF(N7-N7/2*2) 23,24,23
   23 IF(IBUG2-1) 11,25,25
   25 WRITE(IWRITE,104)
      GO TO 11
   24 N1=2*IHSH-1
      N2=IHSH+1
      N3=IHSH-1
      N4=IHSH-2
      GO TO 3
*
* --- THE TWO CONFIGURATIONS WILL HAVE ZERO HAMILTONIAN MATRIX ELEMENT
*
   11 LET=0
      RETURN
    3 CONTINUE
*
* --- NO OBVIOUS ANGULAR MOMENTUM ORTHOGONALITY
*
   12 LET=1
      IF (IELST.EQ.0.AND. .NOT.INCL) LET = 0
      RETURN
      END
*
*     ------------------------------------------------------------------
*       S H E L L S
*     ------------------------------------------------------------------
*
      SUBROUTINE SHELLS(JA,JB,LET)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH(16,2),J1QN(31,3,2),
     :IJFUL(16)
      COMMON/OCCUPATION/NCG(16,2),ICG(16,2)
      COMMON/JB/IHSHJB,NJB(16),LJB(16),NOSHJB(16),IJFULJB(16)
C      WRITE(6,'(/A)') '  **************  iena  ******** '
C      WRITE(6,'(/A)') '  JA JB '
C      WRITE(6,'(2I5)') JA,JB
      LET=1
      IF(JA.EQ.JB) THEN
      CALL SHELLSAME(JA)
      RETURN
      ENDIF
      IA=NOCCSH(JA)
      IB=NOCCSH(JB)
      IF(IABS(IA-IB).GT.4) THEN
      LET=0
      RETURN
      ENDIF
      IAM=0
      IBM=0
      IGA=0
      IGGA=0
      IHSH=0
      JG1=1
      JG2=1
    1 J1=NOCORB(JG1,JA)
      J2=IJFULJB(JG2)
      IHSH=IHSH+1
C      WRITE(6,'(/A)') '  IHSH '
C      WRITE(6,'(I3)')  IHSH 
C      WRITE(6,'(A)') '  JG1, JG2 '
C      WRITE(6,'(2I5)') JG1,JG2
C      WRITE(6,'(A)') '        J1, J2 '
C      WRITE(6,'(5X,2I5)') J1,J2
      IF(J1.EQ.J2) THEN
        NA=2
        NB=2
        JJ=J1
      ELSEIF(J1.GT.J2) THEN
        IF(IBM.EQ.0) THEN
          NA=1
          NB=2
          JJ=J2
	ELSE
          NA=2
          NB=1
          JJ=J1
	ENDIF
      ELSEIF(J1.LT.J2) THEN
        IF(IAM.EQ.0) THEN
          NA=2
          NB=1
          JJ=J1
	ELSE
          NA=1
          NB=2
          JJ=J2
	ENDIF
      ENDIF
      NCG(IHSH,1)=NA
      NCG(IHSH,2)=NB
      ICG(IHSH,1)=JG1
      ICG(IHSH,2)=JG2
      IF(NA.EQ.1) THEN
        NOSH(IHSH,1)=0
      ELSE
        NOSH(IHSH,1)=NELCSH(JG1,JA)
      ENDIF
      IF(NB.EQ.1) THEN
        NOSH(IHSH,2)=0
      ELSE
        NOSH(IHSH,2)=NOSHJB(JG2)
      ENDIF
C      WRITE(6,'(A)') '  nosh(1),nosh(2) '
C      WRITE(6,'(2I5)') NOSH(IHSH,1),NOSH(IHSH,2)
      NJ(IHSH)=NJCOMP(JJ)
C      WRITE(6,'(A)') '  NJ '
C      WRITE(6,'(I5)') NJ(IHSH)
      LJ(IHSH)=LJCOMP(JJ)
C      WRITE(6,'(A)') '  LJ '
C      WRITE(6,'(I5)') LJ(IHSH)
      IJFUL(IHSH)=JJ
      IF(JG1.LT.IA.OR.JG2.LT.IB) THEN
        IF(JG1.LT.IA) THEN
	  IAM=0
          IF(NA.EQ.2) JG1=JG1+1
        ELSE
          IF(NA.EQ.2)  THEN
	  IF(IGGA.EQ.0) THEN	
	  IAM=1
	  IGG1=IGG1+1
        ENDIF
        ENDIF
        ENDIF
        IF(JG2.LT.IB) THEN
	  IBM=0
          IF(NB.EQ.2) JG2=JG2+1
        ELSE
          IF(NB.EQ.2)  THEN
	  IF(IGGA.EQ.0) THEN	
	  IBM=1
	  IGG1=IGG1+1
        ENDIF
        ENDIF
        ENDIF
        GO TO 1
      ENDIF
      IF(NA.EQ.1) THEN
	 IF(J1.LT.J2)RETURN
	 IF(IAM.EQ.1) RETURN
         IAM=1
	 IBM=1
         IGA=IGA+1
         IF(IGA.EQ.1)GO TO 1
      ENDIF
      IF(NB.EQ.1) THEN
	 IF(J2.LT.J1)RETURN
	 IF(IBM.EQ.1) RETURN
         IAM=1
	 IBM=1
         IGA=IGA+1
         IF(IGA.EQ.1)GO TO 1
      ENDIF
      RETURN
      END
*
*     ------------------------------------------------------------------
*       S H E L L S J B
*     ------------------------------------------------------------------
*
      SUBROUTINE SHELLSJB(JB)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      COMMON/OCCUPATION/NCG(16,2),ICG(16,2)
      COMMON/JB/IHSHJB,NJB(16),LJB(16),NOSHJB(16),IJFULJB(16)
      IB=NOCCSH(JB)
      IHSHJB=IB
      DO 1 J=1,IB
        I1=NOCORB(J,JB)
        NOSHJB(J)=NELCSH(J,JB)
        NJB(J)=NJCOMP(I1)
        LJB(J)=LJCOMP(I1)
        IJFULJB(J)=I1
    1 CONTINUE
      RETURN
      END
*
*     ------------------------------------------------------------------
*       S H E L L S A M E
*     ------------------------------------------------------------------
*
      SUBROUTINE SHELLSAME(JA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/JB/IHSHJB,NJB(16),LJB(16),NOSHJB(16),IJFULJB(16)
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH(16,2),J1QN(31,3,2),
     :IJFUL(16)
      COMMON/OCCUPATION/NCG(16,2),ICG(16,2)
      DO 2 J=1,IHSHJB
        NCG(J,1)=2
        NCG(J,2)=2
        ICG(J,1)=J
        ICG(J,2)=J
        NOSH(J,1)=NOSHJB(J)
        NOSH(J,2)=NOSHJB(J)
        NJ(J)=NJB(J)
        LJ(J)=LJB(J)
        IJFUL(J)=IJFULJB(J)
    2 CONTINUE
      IHSH=IHSHJB
      RETURN
      END
*
*     ------------------------------------------------------------------
*      C O U P L I N G 
*     ------------------------------------------------------------------
*
      SUBROUTINE COUPLING(JA,JB)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH(16,2),J1QN(31,3,2),
     :IJFUL(16)
      COMMON/OCCUPATION/NCG(16,2),ICG(16,2)
      DO 1 IH=1,IHSH
        JC=JA
        DO 2 I=1,2
          IC=ICG(IH,I)
          NC=NCG(IH,I)
          I2H=IHSH+IH-1
*
* --- FIRST CONSIDER THE L.H.S. (I=1) OF THE MATRIX ELEMENT. NC=1 MEANS
*     UNOCCUPIED, REPRESENTED BY A DUMMY SINGLET S SHELL, AND THE
*    ADDITIONAL SET OF COUPLING QUANTUM NUMBERS WILL BE THE SAME AS THE
*     LAST SET OF COUPLING QUANTUM NUMBERS ALREADY OBTAINED.
*     NC=2 MEANS OCCUPIED.  THEN ALL THE NEW QUANTUM NUMBERS (BOTH FOR
*     THE SHELL AND FOR THE COUPLING OF THIS SHELL TO THE RESULTANT OF
*     THE PREVIOUS ONES) ARE DEFINED IN THE CORRESPONDING J1QNRD ARRAY.
*     NOSH - THE NUMBER OF ELECTRONS IN THIS SHELL, IS DEFINED BY THE
*     APPROPRIATE ENTRY IN NELCSH .  THE R.H.S. IS THEN CONSIDERED
*     SIMILARLY (I=2)
*
          IF(NC.EQ.1) THEN
            J1QN(IH,1,I)=0
            J1QN(IH,2,I)=1
            J1QN(IH,3,I)=1
            IF(IH.EQ.2) THEN
              J1QN(I2H,1,I)=0
              J1QN(I2H,2,I)=J1QN(1,2,I)
              J1QN(I2H,3,I)=J1QN(1,3,I)
            ELSEIF(IH.GT.2) THEN
	      I2H1=I2H-1
              J1QN(I2H,1,I)=J1QN(I2H1,1,I)
              J1QN(I2H,2,I)=J1QN(I2H1,2,I)
              J1QN(I2H,3,I)=J1QN(I2H1,3,I)
            END IF
          ELSE
          JD = J1QNRD(IC,JC)
          J1QN(IH,1,I)=MOD(JD,64)
          JD = JD/64
          J1QN(IH,2,I) = MOD(JD,64)
          J1QN(IH,3,I) = JD/64
*
*     IS THIS THE FIRST OCCUPIED SHELL OF EITHER CONFIGURATION. IF SO,
*    THEN THERE ARE NO INTERMEDIATE COUPLINGS TO CONSIDER AT THIS STAGE
*
          IF(IH .GT. 1) THEN
*
*    IS THIS THE FIRST OCCUPIED SHELL OF THIS CONFIGURATION, THOUGH NOT
*     THE FIRST OF THE OTHER CONFIGURATION.  IF SO, THE INTERMEDIATE
*     COUPLING FORMED HAS THE SAME  L,S  VALUES AS THIS OCCUPIED SHELL,
*     SINCE WE COUPLE THE SHELL TO A DUMMY SINGLET S.
*
            IF(IC .LE.1) THEN
              I2 = 1
            ELSE
              I2 = NOCCSH(JC)+IC-1
            END IF
            JD = J1QNRD(I2,JC)
            IF (IC .LE. 1) THEN
              J1QN(I2H,1,I) = 0
            ELSE
              J1QN(I2H,1,I) = MOD(JD,64)
            END IF
            JD = JD/64
            J1QN(I2H,2,I) = MOD(JD, 64)
            J1QN(I2H,3,I) = JD/64
          END IF
        END IF
        JC=JB
    2   CONTINUE
    1 CONTINUE
      END
*
*     -------------------------------------------------------------
*      S A V E N O N
*     -------------------------------------------------------------
*                                                                  *
*     THIS SUBROUTINE FOR         G I S O                          *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*
      SUBROUTINE SAVENON(I,A,KL,LA,LB,LC,LD,JA,JB,IPTR)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      PARAMETER (NWD=60)
      COMMON/ANGCORE/CBANG(NWD,NWD)
      COMMON/GRADINT/GGRAD(NWD,NWD)
      POINTER (QWT,WT(1))
      COMMON /NAN/ QWT
      COMMON/ISO/GIS,RIS,COR,CO,FIELCOR,FIEL
      COMMON/SIGNF/SIGNFA
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/CLOSED/B1ELC(4),NCLOSD,IBK
      ILB=LB+NCLOSD
      ILD=LD+NCLOSD
      A=A*SIGNFA
      IF(I.NE.4) THEN
        IF(KL.EQ.1) THEN
          IF(I.EQ.2) THEN
            CISO=WT(JA)*WT(JB)*A*GGRAD(ILB,ILD)**2
            IF(JA.NE.JB) CISO=TWO*CISO
            GIS=GIS+CISO
          ELSEIF(I.EQ.3) THEN
            ILA=LA+NCLOSD
            ILC=LC+NCLOSD
            CISO=-WT(JA)*WT(JB)*A*GGRAD(ILA,ILC)*GGRAD(ILB,ILD)
            IF(JA.NE.JB) CISO=TWO*CISO
            RIS=RIS+CISO
          ENDIF
        ENDIF
      ELSE
        CISO=WT(JA)*WT(JB)*A
        IF(JA.NE.JB) CISO=TWO*CISO
	IF(LB.LE.LD) THEN
          CBANG(LB,LD)=CISO+CBANG(LB,LD)
        ELSE
          CBANG(LD,LB)=CISO+CBANG(LD,LB)
        ENDIF
      ENDIF
      RETURN
      END
