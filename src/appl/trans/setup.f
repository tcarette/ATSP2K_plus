*
*     ------------------------------------------------------------------
*	S E T U P
*     ------------------------------------------------------------------
*
*	In order to keep NC in COMMON /RED/ as in other
*   routines, the role of NCC and NC has been reversed from
*   the earlier version of this code.
*
      SUBROUTINE setup(JA,JB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc(2),iuw(2),iul(2),iuj(2)
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,MAXORB
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH(16,2),J1QN(31,3,2),IJFUL(16)
      POINTER (QIORTH,IORTH(1))
      COMMON/OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
     : ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP,
     : QIORTH
      COMMON/RED/INDL(20),INDR(20),IORST(20),NCI1,KPL,KPR,NC,LIS(16),
     : JIST,JFST,NCONTR
*
*     NOTICE THE DIFFERENT NAMES IN THE COMMON BLOCK MEDEFN  -  WE
*      STORE NOSH1(I=1,10) AS NOSH((I=1,10),1) AND NOSH2(I=1,10) AS
*     NOSH((I=1,10),2)   AND USE THE FACT THAT NOSH1 AND NOSH2 WILL THEN
*     BE EQUIVALENT TO THE SINGLE 2-DIMENSIONAL ARRAY NOSH.  SIMILARLY
*     FOR J1QN
*
* === GENERATES THE ARRAYS  NJ,LJ - DEFINING THE QUANTUM NUMBERS OF THE
*     SHELLS,   NOSH - DEFINING THE OCCUPATION OF THE SHELLS,  J1QN -
*     DEFINING THE COUPLING OF THE SHELLS,   FOR EACH OF THE TWO
*     CONFIGURATIONS CONSIDERED.    ONLY THOSE SHELLS OCCURRING IN AT
*     LEAST ONE CONFIGURATION ARE INCLUDED.
*                   AT LEAST TWO SHELLS MUST BE CONSIDERED OCCUPIED.
*     THUS (1S)**2    HELIUM  MUST BE TREATED AS ,E.G., (1S)**2(2S)**0
*     THE SIZE OF THE ARRAYS HERE CALCULATED IS ARRANGED TO BE NO
*     GREATER THAN IS NECESSARY TO INCLUDE ALL ORBITALS WHICH ARE
*     DEEMED TO BE OCCUPIED IN EITHER OR BOTH OF THE CONFIGURATIONS
*     JA,JB
*
* --- INITIALIZE BASIC QUANTITIES - (I1+1) RUNS OVER 1,MAXORB,  IHSH IS
*     THE CURRENT VALUE OF THE HIGHEST OCCUPIED SHELL YET CONSIDERED,
*     WHILE I2HSH=2*IHSH-1
*
      MU=0
      NU=0
      MUP=0
      NUP=0
      I1=0
      IHSH=0
      I2HSH=-1
      IA=NOCCSH(JA)
      IB=NOCCSH(JB)
*
* --- TEST ON WHETHER LIMIT OF I1 HAS BEEN REACHED
*
    1 IF(I1-MAXORB) 101,100,100
*
* --- INCREASE BASIC QUANTITIES
*
  101 I1=I1+1
      I3=IHSH+1
      I5=I2HSH+I3
*
* --- IS THE SHELL I1 OCCUPIED IN JA
*
      DO 2 J=1,IA
      IF(I1-NOCORB(J,JA)) 2,3,2
    2 CONTINUE
      NA=1
      GO TO 4
    3 NA=2
      J1=J
*
* --- IS THE SHELL I1 OCCUPIED IN JB
*
    4 DO 5 J=1,IB
      IF(I1-NOCORB(J,JB)) 5,6,5
    5 CONTINUE
      NB=1
      GO TO 7
    6 NB=2
      J2=J
*
*     IF THE SHELL I1 IS NOT OCCUPIED IN EITHER JA OR JB, IGNORE THE
*     SHELL, DO NOT INCREASE IHSH, AND CONSIDER NEXT SHELL BY INCREASING
*     I1
*
    7 IF(NA-1) 8,8,9
    8 IF(NB-1) 1,1,9
*
* --- IF THE SHELL I1 IS OCCUPIED IN EITHER JA OR JB -
*     (1)   IF IHSH.GT.1, THEN ALREADY AT LEAST TWO SHELLS AND THE
*     RESULTING COUPLINGS HAVE BEEN STORED. WE MUST THUS MAKE ROOM FOR
*     THE QUANTUM NUMBERS OF THIS NEW SHELL BETWEEN THE QUANTUM NUMBERS
*     OF THE PREVIOUS SHELLS AND THE QUANTUM NUMBERS OF THE INTERMEDIATE
*     COUPLINGS OF THE CONFIGURATIONS.  THUS THE LATTER SET ARE =MOVED
*     ALONG= TO MAKE ROOM FOR THE NEW SHELL
*     (2)   IF IHSH.LE.1, THERE ARE NO INTERMEDIATE COUPLING QUANTUM
*     NUMBERS, AND SO THERE IS NOTHING TO MOVE
*
    9 IF(IHSH-1) 11,11,10
   10 DO 12 I=1,2
      DO 13 J=I3,I2HSH
      I4=I5-J
      DO 14 K=1,3
      J1QN(I4+1,K,I)=J1QN(I4,K,I)
   14 CONTINUE
   13 CONTINUE
   12 CONTINUE
   11 IHSH=I3
      I2HSH=I2HSH+2
      NCC=NA
      I=1
      IC=J1
      JC=JA
*
* --- FIRST CONSIDER THE L.H.S. (I=1) OF THE MATRIX ELEMENT.  NC=1 MEANS
*     UNOCCUPIED, REPRESENTED BY A DUMMY SINGLET S SHELL, AND THE
*     ADDITIONAL SET OF COUPLING QUANTUM NUMBERS WILL BE THE SAME AS THE
*     LAST SET OF COUPLING QUANTUM NUMBERS ALREADY OBTAINED.
*     NC=2 MEANS OCCUPIED.  THEN ALL THE NEW QUANTUM NUMBERS (BOTH FOR
*     THE SHELL AND FOR THE COUPLING OF THIS SHELL TO THE RESULTANT OF
*     THE PREVIOUS ONES) ARE DEFINED IN THE CORRESPONDING J1QNRD ARRAY.
*     NOSH - THE NUMBER OF ELECTRONS IN THIS SHELL, IS DEFINED BY THE
*     APPROPRIATE ENTRY IN NELCSH .  THE R.H.S. IS THEN CONSIDERED
*     SIMILARLY (I=2)
*
   25 GO TO (15,16),NCC
   15 NOSH(IHSH,I)=0
      J1QN(IHSH,1,I)=0
      J1QN(IHSH,2,I)=1
      J1QN(IHSH,3,I)=1
      IF(IHSH-2) 22,18,19
   18 J1QN(3,1,I)=0
      J1QN(3,2,I)=J1QN(1,2,I)
      J1QN(3,3,I)=J1QN(1,3,I)
      GO TO 22
   19 DO 27 K=1,3
      J1QN(I2HSH,K,I)=J1QN(I2HSH-1,K,I)
   27 CONTINUE
      GO TO 22
   16 NOSH(IHSH,I)=NELCSH(IC,JC)
      IF(NOVLPS.EQ.0) GO TO 33
      GO TO (31,32),I
   31 IF(IC.EQ.JMU) MU=IHSH
      IF(IC.EQ.JNU) NU=IHSH
      GO TO 33
   32 IF(IC.EQ.JMUP) MUP=IHSH
      IF(IC.EQ.JNUP) NUP=IHSH
   33 JD = J1QNRD(IC,JC)
      J1QN(IHSH,1,I)=MOD(JD,64)
      JD = JD/64
      J1QN(IHSH,2,I) = MOD(JD,64)
      J1QN(IHSH,3,I) = JD/64
*
*     IS THIS THE FIRST OCCUPIED SHELL OF EITHER CONFIGURATION. IF SO,
*     THEN THERE ARE NO INTERMEDIATE COUPLINGS TO CONSIDER AT THIS STAGE
*
      IF(IHSH .GT. 1) THEN
*
*     IS THIS THE FIRST OCCUPIED SHELL OF THIS CONFIGURATION, THOUGH NOT
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
	    J1QN(I2HSH,1,I) = 0
         ELSE
	    J1QN(I2HSH,1,I) = MOD(JD,64)
         END IF
         JD = JD/64
         J1QN(I2HSH,2,I) = MOD(JD, 64)
         J1QN(I2HSH,3,I) = JD/64
      END IF
*
*     SENIORITY SET (ARBITRARILY) ZERO FOR INTERMEDIATE COUPLING
*
   22 IF(I-2) 23,24,24
   23 NCC=NB
      I=2
      IC=J2
      JC=JB
      GO TO 25
*
* --- SET THE NJ AND LJ VALUES OF THE OCCUPIED SHELLS
*
   24 NJ(IHSH)=NJCOMP(I1)
      IJFUL(IHSH)=I1
      LJ(IHSH)=LJCOMP(I1)
      IF (NC .NE. 0) THEN
         IF (INDL(NCI1).EQ.I1) KPL=IHSH
         IF (INDR(NCI1).EQ.I1) KPR=IHSH
      ENDIF
*
* --- RETURN TO 1  TO SEE IF MAXORB HAS BEEN REACHED
*
      GO TO 1
  100 RETURN
      END
