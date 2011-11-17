*
*     -------------------------------------------------------------
*      T W O P A R T I C L E 3 3
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
*                                              N'2 = N2 + 1        *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE TWOPARTICLE33(IA,IB,IC,IIA,IIB,IIC,IID)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      COMMON/TRK2/BD3(3),BD4(3),BK3(3),BK4(3),
     *ID3(7),ID4(7),IK3(7),IK4(7)
      COMMON/PERMAT/IRS(2,2),IRL(20,20),RS(2,2),RL(20,20)
      COMMON/KAMPAS/ IW1(2,20),IW2(2,20),IWAA(2,2,20,20),
     :RW1(2,20),RW2(2,20),RWAA(2,2,20,20)
      COMMON /CASEOP/ IOCASE
      EXTERNAL TWO33A,TWO33B
      IF(IHSH.LE.2)RETURN
      CALL HIBFF(IA,IB,IC,IA,3)
      IF(ITTK(1,IK1(6),ID1(6)).EQ.0) RETURN
      IF(ITTK(1,IK2(6),ID2(6)).EQ.0) RETURN
      IF(IABS(IK3(6)-ID3(6)).GT.2) RETURN
      IF(ITTK(2*IK1(3),IK1(5),ID1(5)).EQ.0) RETURN
      IF(ITTK(2*IK2(3),IK2(5),ID2(5)).EQ.0) RETURN
      IF(IABS(IK3(5)-ID3(5)).GT.4*IK3(3)) RETURN
      DO 3 I1=1,2
        DO 4 I2=1,2
	  IRS(I1,I2)=0
    4   CONTINUE
    3 CONTINUE
      DO 5 I1=1,20
	IW1(1,I1)=0
	IW1(2,I1)=0
        DO 6 I2=1,20
	  IRL(I1,I2)=0
    6   CONTINUE
    5 CONTINUE
      IF(IOCASE.EQ.1) THEN
* Spin-Spin operator
        IP=ITREXG2(LJ(IA),LJ(IB),LJ(IC),LJ(IC),IKK)+1
      ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
        IP=ITREXG(LJ(IA),LJ(IB),LJ(IC),LJ(IC),IKK)+1
      ENDIF
      IF(IKK.GT.0) THEN
        IG=IP+IKK-1
        DO 1 I=IP,IG
          KL=I-1
*
          IF(IOCASE.EQ.1) THEN
* Spin-spin operator
            CALL SSC(IG,KL,IA,IB,IC,IC,IB,IC,IA,IC,KL,TWO33A)
          ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
            CALL SOOC(IG,KL,IA,IB,IC,IC,IB,IC,IA,IC,KL,TWO33A)
          ENDIF
*
    1   CONTINUE
      ENDIF
      IF(IOCASE.EQ.1) THEN
* Spin-Spin operator
        IP=ITREXG2(LJ(IA),LJ(IC),LJ(IB),LJ(IC),IKK)+1
      ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
        IP=ITREXG(LJ(IA),LJ(IC),LJ(IB),LJ(IC),IKK)+1
      ENDIF
      IF(IKK.LE.0)RETURN
      IG=IP+IKK-1
      DO 2 I=IP,IG
        KL=I-1
        IF(IOCASE.EQ.1) THEN
* Spin-spin operator
          CALL SSC(IG,KL,IA,IB,IC,IC,IC,IB,IA,IC,KL,TWO33B)
        ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
          CALL SOOC(IG,KL,IA,IB,IC,IC,IC,IB,IA,IC,KL,TWO33B)
	ENDIF
    2 CONTINUE
      RETURN
      END
