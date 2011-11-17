*
*     -------------------------------------------------------------
*      T W O 2 A 
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
      SUBROUTINE TWO2A(IA,IB)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      COMMON/PERMAT/IRS(2,2),IRL(20,20),RS(2,2),RL(20,20)
      COMMON/KAMPAS/ IW1(2,20),IW2(2,20),IWAA(2,2,20,20),
     :RW1(2,20),RW2(2,20),RWAA(2,2,20,20)
      COMMON /CASEOP/ IOCASE
      EXTERNAL TWO2
      IF(IHSH.LE.1)RETURN
*
      IF(IOCASE.EQ.1) THEN
* Spin - spin      operator
	KGAL=2
      ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
	KGAL=1
      ENDIF
*
      IIA=MIN0(IA,IB)
      IIB=MAX0(IA,IB)
      CALL RLSP0(1,IIA,IIB,0,IAT)
      IF(IAT.EQ.0) RETURN
      CALL RLSP0(2,IIA,IIB,2*KGAL,IAT)
      IF(IAT.EQ.0) RETURN
      CALL RLSP0(3,IIA,IIB,2*KGAL,IAT)
      IF(IAT.EQ.0) RETURN
      CALL HIBFF(IA,IB,IA,IA,2)
      IF(IABS(IK1(6)-ID1(6)).GT.2) RETURN
      IF(IABS(IK2(6)-ID2(6)).GT.2) RETURN
      IF(IABS(IK1(5)-ID1(5)).GT.4*IK1(3)) RETURN
      IF(IABS(IK2(5)-ID2(5)).GT.4*IK2(3)) RETURN
      DO 3 I1=1,2
        DO 4 I2=1,2
          IRS(I1,I2)=0
    4   CONTINUE
    3 CONTINUE
      DO 5 I1=1,20
        IW1(1,I1)=0
        IW1(2,I1)=0
        IW2(1,I1)=0
        IW2(2,I1)=0
        DO 6 I2=1,20
          IRL(I1,I2)=0
    6   CONTINUE
    5 CONTINUE
      IP=IABS(LJ(IB)-LJ(IA))+1
      IG=LJ(IB)+LJ(IA)+1
      DO 1 I=IP,IG,2
        KL=I-1
        IF(IOCASE.EQ.1) THEN
* Spin-spin operator
          CALL SS1122(IG,KL,IA,IB)
        ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
	  CALL SOO1122P(IG,KL,IA,IB,IB,IB,IA,IA,IB,IB,KL,TWO2)
C          CALL SOO1122(IG,KL,IA,IB,TWO2)
        ENDIF
    1 CONTINUE
      RETURN
      END
