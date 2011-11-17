*                                                                  : 
*     ..........................................................   :
*        ii)   two - particle operator                             : 
*     ..........................................................   :
*
*     -------------------------------------------------------------
*      T W O 1 
*     -------------------------------------------------------------
*                                                                  *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
*                                                                  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE TWO1(IIA,IIB,IIRE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      COMMON/PERMAT/IRS(2,2),IRL(20,20),RS(2,2),RL(20,20)
      COMMON/KAMPAS/ IW1(2,20),IW2(2,20),IWAA(2,2,20,20),
     :RW1(2,20),RW2(2,20),RWAA(2,2,20,20)
      COMMON /CASEOP/ IOCASE
      EXTERNAL TWO13
*
      IF(IOCASE.EQ.1) THEN
* Spin - spin      operator
	KGAL=2
      ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
	KGAL=1
      ENDIF
*
      IF(IIA.EQ.IIB) THEN
        IA=IIA
        IF(ITTK(J1QN1(IA,2)-1,J1QN2(IA,2)-1,2*KGAL).EQ.0) RETURN
        IF(ITTK(J1QN1(IA,3)-1,J1QN2(IA,3)-1,2*KGAL).EQ.0) RETURN
        CALL HIBFF(IA,IA,IA,IA,1)
        IF(ID1(4).LT.2)RETURN
        CALL RLSP0(1,IA,IA,0,IAT)
        IF(IAT.EQ.0) RETURN
        CALL RLSP0(2,IA,IA,2*KGAL,IAT)
        IF(IAT.EQ.0) RETURN
        CALL RLSP0(3,IA,IA,2*KGAL,IAT)
        IF(IAT.EQ.0) RETURN
        CALL RLSP1(3,IA,2*KGAL,0,IAT,REC)
        IF(IAT.EQ.0) RETURN
        CALL RLSP1(2,IA,2*KGAL,0,IAT,REC)
        IF(IAT.EQ.0) RETURN
        CALL RLSP1(3,IA,2*KGAL,1,IAT,RECS)
        IF(DABS(RECS).LT.EPS) RETURN
        CALL RLSP1(2,IA,2*KGAL,1,IAT,RECL)
        IF(DABS(RECL).LT.EPS) RETURN
	RS(1,1)=RECS
	RL(1,1)=RECL
        LIA2=LJ(IA)*2
        IG=LIA2+1
	IRS(1,1)=0
        DO 1 I=1,IG,2
          KL=I-1
          IF(IOCASE.EQ.1) THEN
* Spin - spin      operator
            CALL SS1111(IG,KL,IA)
          ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
            CALL SOO1111(IG,KL,IA)
C            CALL SOO1111P(IG,KL,IA)
          ENDIF
    1   CONTINUE
      ELSE
        IF(IIRE.EQ.0)RETURN
        IF(IHSH.LE.1)RETURN
        IA=MIN0(IIA,IIB)
        IB=MAX0(IIA,IIB)
        CALL RLSP0(1,IA,IB,0,IAT)
        IF(IAT.EQ.0) RETURN
        CALL RLSP0(2,IA,IB,2*KGAL,IAT)
        IF(IAT.EQ.0) RETURN
        CALL RLSP0(3,IA,IB,2*KGAL,IAT)
        IF(IAT.EQ.0) RETURN
        CALL HIBFF(IA,IB,IA,IA,2)
        IF(IABS(IK1(6)-ID1(6)).GT.2) RETURN
        IF(IABS(IK2(6)-ID2(6)).GT.2) RETURN
        IF(IABS(IK1(5)-ID1(5)).GT.4*IK1(3)) RETURN
        IF(IABS(IK2(5)-ID2(5)).GT.4*IK2(3)) RETURN
        DO 4 I1=1,2
          DO 5 I2=1,2
	    IRS(I1,I2)=0
    5     CONTINUE
    4   CONTINUE
        DO 6 I1=1,20
	  IW1(1,I1)=0
	  IW1(2,I1)=0
	  IW2(1,I1)=0
	  IW2(2,I1)=0
          DO 7 I2=1,20
	    IRL(I1,I2)=0
    7     CONTINUE
    6   CONTINUE
        LIA2=LJ(IA)*2
        LIB2=LJ(IB)*2
        IG=MIN0(LIA2,LIB2)+1
        DO 2 I=1,IG,2
          KL=I-1
          IF(IOCASE.EQ.1) THEN
* Spin - spin      operator
            CALL SS1212(IG,KL,IA,IB)
          ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
            CALL SOO1212(IG,KL,IA,IB)
C            CALL SOO1212P(IG,KL,IA,IB)
          ENDIF
    2   CONTINUE
        IP=IABS(LJ(IB)-LJ(IA))+1
        IG=LJ(IB)+LJ(IA)+1
        DO 3 I=IP,IG,2
          KL=I-1
          IF(IOCASE.EQ.1) THEN
* Spin-spin operator
            CALL SS1221(IG,KL,IA,IB)
          ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
            CALL SOOC(IG,KL,IA,IB,IB,IB,IA,IB,IB,IA,KL,TWO13)
          ENDIF
    3   CONTINUE
      ENDIF
      RETURN
      END
