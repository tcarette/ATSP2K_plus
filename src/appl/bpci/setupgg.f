*     ..........................................................   :
*                                                                  : 
*          Block                                                   : 
*                             S E T U P G G                        :
*                                                                  : 
*                                                                  : 
*     Written by G. Gaigalas,                                      : 
*                  Department  of  Computer Science,               : 
*                  Vanderbilt University,  Nashville               : 
*                                                  February 1994   : 
*                                                                  : 
*     ..........................................................   :
*
*
*     ----------------------------------------------------------------
*      S E T U P G G  
*     ----------------------------------------------------------------
*
*     THE ROUTINE EVALUATES THE NON-RELATIVISTIC HAMILTONIAN
*     WITH ORTHOGONAL ORBITALS
*
*
*     THE MATRIX ELEMENT OF THE TWO-ELECTRON POTENTIAL BETWEEN TWO
*     STATES (LABELLED 1 AND 2) MAY BE EXPRESSED AS A SUM OF WEIGHTED
*     RK (SLATER) INTEGRALS.  THIS SUBROUTINE, TOGETHER WITH THOSE
*     CALLED BY IT, DETERMINES THESE WEIGHTS, WHICH ARISE FROM AN
*     INTEGRATION OVER THE ANGULAR AND SPIN CO-ORDINATES
*     FOR DETAILS, SEE   U. FANO, PHYS. REV.,140,A67,(1965)
*
*     THE =INTERACTING= SHELLS ARE DESIGNATED  IRHO,ISIG,IRHOP,ISIGP.
*     THE FIRST TWO REFER TO THE L.H.S. OF     (PSI/V/PSIP)     , WHILE
*     THE SECOND TWO REFER TO THE R.H.S.  FOR DIAGONAL AND CERTAIN OFF-
*     DIAGONAL MATRIX ELEMENTS, THESE MAY NOT BE UNIQUE, AND EACH
*     POSSIBILITY MUST BE CONSIDERED IN TURN
*     THE CONDITION =IRHO .LE. ISIG ,  IRHOP .LE. ISIGP=  IS TO BE
*     SATISFIED
*
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*
      SUBROUTINE SETUPGG
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(8),ISCW
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/DIAGNL/IDIAG,JA,JB
      COMMON/STEGG/IX,IGGG,IRHO,ISIG,IRHOP,ISIGP
    5 FORMAT(//10X,7H IRHO =,I3,4X,7H ISIG =,I3,4X,8H IRHOP =,I3,3X,8H I
     :SIGP =,I3)
      IZERO = 0
      IX=0
      IRHO=0
      ISIG=0
      IRHOP=0
      ISIGP=0
      IGGG=0
      DO 1 J=1,IHSH
        N=NOSH1(J)-NOSH2(J)
        IF(IABS(N).GT.2) THEN
*         ... added by Anders Ynnerman, IX needs to be set upon return. 
*                                       ..................
*                                       Thank You very much. 
*                                                           Gediminas 
          IX = 5
          RETURN
        END IF
        IF(N.GT.0) THEN
          IF(N.EQ.1) THEN
            ISIG=J
            IF(IRHO.EQ.0) IRHO = J
            IX =IX+1
          ELSE
            IRHO=J
            IF(IABS(N).EQ.2)IGGG=20+IGGG
            IX=IX+2
          ENDIF
        ELSEIF(N.LT.0) THEN
          IF(N+1.EQ.0) THEN
            ISIGP=J
            IF(IRHOP.EQ.0) IRHOP = J
            IX=IX+1
          ELSE
            IRHOP=J
            IF(IABS(N).EQ.2)IGGG=20+IGGG
            IX=IX+2
          ENDIF
        ENDIF
    1 CONTINUE
*
*     IX MEASURES THE TOTAL NUMBER OF ELECTRONS IN EITHER CONFIGURATION
*     WHICH DO NOT OCCUR IN THE OTHER.  THEN IF  IX  IS GREATER THAN 4,
*     ORTHOGONALITY OF THE ORBITALS PREVENTS A NON-ZERO MATRIX ELEMENT.
*     IF  IX  IS LESS THAN 4, THEN WE DIVIDE IX BY 2 AND NOW IX MEASURES
*     THE NUMBER OF ELECTRONS WHICH HAVE BEEN CHANGED IN GOING FROM PSI
*     TO PSIP.  IF NOW IX=0, WE HAVE A DIAGONAL MATRIX ELEMENT.  RHO AND
*     SIG MAY TAKE ON ANY VALUES LESS THAN IHSH.  IF IX=1, ONE INTER-
*     ACTING SHELL ON EACH SIDE IS FIXED, WHILE THE OTHER MAY VARY.  IF
*     IX=2, ALL INTERACTING SHELLS ARE DETERMINED
*
      IF(IX.GT.4) RETURN
      CALL COUPLING(JA,JB)
      RETURN
      END
*
*     ----------------------------------------------------------------
*      N O N B P
*     ----------------------------------------------------------------
*
*     THE ROUTINE EVALUATES THE NON-RELATIVISTIC HAMILTONIAN
*     WITH ORTHOGONAL ORBITALS
*
*
*     THE MATRIX ELEMENT OF THE TWO-ELECTRON POTENTIAL BETWEEN TWO
*     STATES (LABELLED 1 AND 2) MAY BE EXPRESSED AS A SUM OF WEIGHTED
*     RK (SLATER) INTEGRALS.  THIS SUBROUTINE, TOGETHER WITH THOSE
*     CALLED BY IT, DETERMINES THESE WEIGHTS, WHICH ARISE FROM AN
*     INTEGRATION OVER THE ANGULAR AND SPIN CO-ORDINATES
*     FOR DETAILS, SEE   U. FANO, PHYS. REV.,140,A67,(1965)
*
*     THE =INTERACTING= SHELLS ARE DESIGNATED  IRHO,ISIG,IRHOP,ISIGP.
*     THE FIRST TWO REFER TO THE L.H.S. OF     (PSI/V/PSIP)     , WHILE
*     THE SECOND TWO REFER TO THE R.H.S.  FOR DIAGONAL AND CERTAIN OFF-
*     DIAGONAL MATRIX ELEMENTS, THESE MAY NOT BE UNIQUE, AND EACH
*     POSSIBILITY MUST BE CONSIDERED IN TURN
*     THE CONDITION =IRHO .LE. ISIG ,  IRHOP .LE. ISIGP=  IS TO BE
*     SATISFIED
*
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*
      SUBROUTINE NONBP
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(8),ISCW
CGG      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(8)
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      COMMON/TRK2/BD3(3),BD4(3),BK3(3),BK4(3),
     *ID3(7),ID4(7),IK3(7),IK4(7)
      COMMON/DIAGNL/IDIAG,JA,JB
      COMMON /OPERAT/ ICOLOM,ISOTOP,IORBORB
      COMMON /BREIT/ ISPORB,ISOORB,ISPSPN
      COMMON/STEGG/IXX,IGGG,IRHO,ISIG,IRHOP,ISIGP
      COMMON /CASEOP/ IOCASE
      EXTERNAL SPINOR
    5 FORMAT(//10X,7H IRHO =,I3,4X,7H ISIG =,I3,4X,8H IRHOP =,I3,3X,8H I
     :SIGP =,I3)
      IZERO = 0
      N1=2*IHSH-1
      JCOLOM=ICOLOM
      JSO=ISPORB
      JSOO=ISOORB
      JSS=ISPSPN
COO
      JORBORB=IORBORB
      IF(ITTK(J1QN1(N1,2)-1,J1QN2(N1,2)-1,0).EQ.0)JCOLOM=0
      IF(ITTK(J1QN1(N1,3)-1,J1QN2(N1,3)-1,0).EQ.0)JCOLOM=0
      IF(ITTK(J1QN1(N1,2)-1,J1QN2(N1,2)-1,0).EQ.0)JORBORB=0
      IF(ITTK(J1QN1(N1,3)-1,J1QN2(N1,3)-1,0).EQ.0)JORBORB=0
COO
      IF(ITTK(J1QN1(N1,2)-1,J1QN2(N1,2)-1,2).EQ.0) THEN
	JSO=0
        JSOO=0
      ENDIF
      IF(ITTK(J1QN1(N1,3)-1,J1QN2(N1,3)-1,2).EQ.0) THEN
	JSO=0
        JSOO=0
      ENDIF
      IF(ITTK(J1QN1(N1,2)-1,J1QN2(N1,2)-1,4).EQ.0)JSS=0
      IF(ITTK(J1QN1(N1,3)-1,J1QN2(N1,3)-1,4).EQ.0)JSS=0
      IX=IXX/2
      IF(IX-1.GT.0) THEN
*
* === UNIQUE SPECIFICATION OF INTERACTING SHELLS
*
        IF(ISIG.EQ.0) ISIG=IRHO
        IF(ISIGP.EQ.0) ISIGP = IRHOP
        IF(IBUG2.GT.0) WRITE(IWRITE,5) IRHO,ISIG,IRHOP,ISIGP
        IF(IGGG.EQ.40) THEN
          IF(JCOLOM.EQ.1) CALL NONRELAT2(IRHO,IRHOP)
          IF(JSOO.EQ.1) THEN
            IOCASE=2
            CALL TWO2A(IRHO,IRHOP)
          ENDIF
          IF(JSS.EQ.1) THEN
            IOCASE=1
            CALL TWO2A(IRHO,IRHOP)
          ENDIF
        ELSEIF(IGGG.EQ.20) THEN
          IF(JCOLOM.EQ.1) CALL NONRELAT4(IRHO,ISIG,IRHOP,ISIGP)
          IF(JSOO.EQ.1) THEN
            IOCASE=2
            CALL TWOPARTICLE4(IRHO,ISIG,IRHOP,ISIGP)
          ENDIF
          IF(JSS.EQ.1) THEN
            IOCASE=1
            CALL TWOPARTICLE4(IRHO,ISIG,IRHOP,ISIGP)
          ENDIF
        ELSE
C          CALL HIBFF(IA,IB,IC,ID,4)
          CALL HIBFF(IRHO,ISIG,IRHOP,ISIGP,4)
          IF(ITTK(1,IK1(6),ID1(6)).EQ.0) RETURN
          IF(ITTK(1,IK2(6),ID2(6)).EQ.0) RETURN
          IF(ITTK(1,IK3(6),ID3(6)).EQ.0) RETURN
          IF(ITTK(1,IK4(6),ID4(6)).EQ.0) RETURN
          IF(ITTK(2*IK1(3),IK1(5),ID1(5)).EQ.0) RETURN
          IF(ITTK(2*IK2(3),IK2(5),ID2(5)).EQ.0) RETURN
          IF(ITTK(2*IK3(3),IK3(5),ID3(5)).EQ.0) RETURN
          IF(ITTK(2*IK4(3),IK4(5),ID4(5)).EQ.0) RETURN
          IF(JCOLOM.EQ.1) CALL NONRELAT5(IRHO,ISIG,IRHOP,ISIGP)
          IF(JSOO.EQ.1) THEN
            IOCASE=2
            CALL TWOPARTICLE5(IRHO,ISIG,IRHOP,ISIGP)
          ENDIF
          IF(JSS.EQ.1) THEN
            IOCASE=1
            CALL TWOPARTICLE5(IRHO,ISIG,IRHOP,ISIGP)
          ENDIF
        ENDIF
      ELSEIF(IX-1.EQ.0) THEN
*
* === ONE INTERACTING SHELL SPECIFIED ON EACH SIDE. SUMMATION OVER OTHER
*
        IRSTO=IRHO
        IRPSTO=IRHOP
        IF(JSO.EQ.1) CALL ONEPARTICLE2(1,1,IRHO,IRHOP,SPINOR)
        DO 2 K1=1,IHSH
          IINE=1
          IIRE=1
          IF(NOSH1(K1).EQ.0) IIRE=0
          ISIG=K1
          IF(NOSH2(K1).EQ.0) IIRE=0
          ISIGP = K1
          IRHO=IRSTO
          IRHOP=IRPSTO
*
*     ORTHOGONALITY OF THE ORBITALS REQUIRES THAT THE VARYING INTER-
*     ACTING SHELL BE THE SAME ON BOTH SIDES OF THE MATRIX ELEMENT
*
* --- IRHO.LE.ISIG,   IRHOP.LE.ISIGP
*
          IF(IRHO.GT.ISIG) THEN
            ISTO=IRHO
            IRHO=ISIG
            ISIG = ISTO
          ELSEIF(IRHO.EQ.ISIG) THEN
            IF(NOSH1(ISIG).EQ.1) IIRE=0
            IF(NOSH1(ISIG).EQ.0) IINE=0
          ENDIF
          IF(IRHOP.GT.ISIGP) THEN
            ISTO=IRHOP
            IRHOP = ISIGP
            ISIGP = ISTO
          ELSEIF(IRHOP.EQ.ISIGP) THEN
            IF(NOSH2(ISIGP).EQ.1) IIRE=0
            IF(NOSH2(ISIGP).EQ.0) IINE=0
          ENDIF
          IF(IBUG2.GT.0) WRITE(IWRITE,5) IRHO,ISIG,IRHOP,ISIGP
          IF(IINE.EQ.1) THEN
            IF(JCOLOM.EQ.1)CALL NONRELAT3(IRHO,ISIG,IRHOP,ISIGP,IIRE)
          ENDIF
          IF(IIRE.EQ.1) THEN
            IF(JSOO.EQ.1) THEN
              IOCASE=2
              CALL TWOPARTICLE3(IRHO,ISIG,IRHOP,ISIGP)
            ENDIF
            IF(JSS.EQ.1) THEN
              IOCASE=1
              CALL TWOPARTICLE3(IRHO,ISIG,IRHOP,ISIGP)
            ENDIF
          ENDIF
    2   CONTINUE
      ELSEIF(IX-1.LT.0) THEN
*
* === NO INTERACTING SHELLS SPECIFIED
*     SUMMATION OVER ALL POSSIBLE COMBINATIONS
*     IN THIS CASE, ORTHOGONALITY OF ORBITALS PRECLUDES ALL CASES
*     EXCEPT  IRHO=IRHOP    AND    ISIG=ISIGP
*
        DO 3 K1=1,IHSH
          IF(NOSH1(K1).NE.0) THEN
            ISIG=K1
            DO 4 K2=1,K1
              IIRE=1
              IF(NOSH1(K2).NE.0) THEN
                IRHO=K2
                IF(IRHO.EQ.ISIG) THEN
                  IF(NOSH1(ISIG).EQ.1) IIRE=0
                ENDIF
                IRHOP=IRHO
                ISIGP=ISIG
                IF(IBUG2.GT.0) WRITE(IWRITE,5) IRHO,ISIG,IRHOP,ISIGP
COO
               IF((JCOLOM+JORBORB).NE.0) CALL NONRELAT1(IRHO,ISIG,IIRE)
COO
                IF(JSOO.EQ.1) THEN
                  IOCASE=2
	          CALL TWO1(IRHO,ISIG,1)
                ENDIF
                IF(JSS.EQ.1) THEN
                  IOCASE=1
	          CALL TWO1(IRHO,ISIG,1)
                ENDIF
              ENDIF
    4       CONTINUE
            IF(JSO.EQ.1) CALL ONEPARTICLE1(1,1,K1,SPINOR)
          ENDIF
    3   CONTINUE
      ENDIF
      RETURN
      END
*
*     ----------------------------------------------------------------
*       S H E L L S
*     ----------------------------------------------------------------
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
      NJ(IHSH)=NJCOMP(JJ)
      LJ(IHSH)=LJCOMP(JJ)
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
*     ----------------------------------------------------------------
*       S H E L L S J B
*     ----------------------------------------------------------------
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
*     ----------------------------------------------------------------
*       S H E L L S A M E
*     ----------------------------------------------------------------
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
*     ----------------------------------------------------------------
*      C O U P L I N G 
*     ----------------------------------------------------------------
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
*     THIS SUBROUTINE FOR      G B R E I T,   G B R C I            *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*
      SUBROUTINE SAVENON(I,A,KL,LA,LB,LC,LD,JA,JB,IPTR)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      POINTER (QSIGNFA,SIGNFA(1))
      COMMON/PHASES/QSIGNFA,ICSTAS
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      SIGNCH=ONE
      IF(ICSTAS.NE.0) SIGNCH=SIGNFA(JA)*SIGNFA(JB)
      A=A*SIGNCH
      CALL SAVE(I,A,KL,LA,LB,LC,LD,JA,JB,IPTR)
      RETURN
      END
