*
*     -------------------------------------------------------------
*      T W O P A R T I C L E 4 
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :      N'1 = N1 +- 1        *
*                                             N'2 = N2 +- 1        *
*                                             N'3 = N3 -+ 2        *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*                                                                  *
      SUBROUTINE TWOPARTICLE4(IA,IB,IC,ID)
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
      EXTERNAL TWO41,TWO42
      IF(IHSH.LE.2)RETURN
*
      IF(IOCASE.EQ.1) THEN
* Spin - spin      operator
	KGAL=2
      ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
	KGAL=1
      ENDIF
*
      IF(IA.EQ.IB) THEN
        CALL EILE(IA,IC,ID,IAA,IBB,ICC)
        CALL HIBFF(IC,ID,IA,IC,3)
      ELSEIF(IC.EQ.ID) THEN
        CALL EILE(IA,IB,IC,IAA,IBB,ICC)
        CALL HIBFF(IA,IB,IC,IA,3)
        IF(ID3(4).LT.2)RETURN
      ELSE
        WRITE(6,'(A)')  ' ERROR IN SUBROUTINE TWOPARTICLE4  '
        STOP
      ENDIF
      IF(ITTK(1,IK1(6),ID1(6)).EQ.0) RETURN
      IF(ITTK(1,IK2(6),ID2(6)).EQ.0) RETURN
      IF(IABS(IK3(6)-ID3(6)).GT.2) RETURN
      IF(ITTK(2*IK1(3),IK1(5),ID1(5)).EQ.0) RETURN
      IF(ITTK(2*IK2(3),IK2(5),ID2(5)).EQ.0) RETURN
      IF(IABS(IK3(5)-ID3(5)).GT.4*IK3(3)) RETURN
      CALL RLSP00(1,IAA,IBB,ICC,ICC,0,IAT)
      IF(IAT.EQ.0) RETURN
      CALL RLSP00(2,IAA,IBB,ICC,ICC,2*KGAL,IAT)
      IF(IAT.EQ.0) RETURN
      CALL RLSP00(3,IAA,IBB,ICC,ICC,2*KGAL,IAT)
      IF(IAT.EQ.0) RETURN
      IF(IOCASE.EQ.1) THEN
* Spin-Spin operator
        IP=ITREXG2(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK)+1
C        IP=ITREXG2(LJ(IA),LJ(IC),LJ(IB),LJ(IC),IKK)+1
      ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
        IP=ITREXG(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK)+1
C        IP=ITREXG(LJ(IA),LJ(IC),LJ(IB),LJ(IC),IKK)+1
      ENDIF
      IF(IKK.LE.0)RETURN
      IG=IP+IKK-1
      DO 2 I1=1,2
        DO 3 I2=1,2
          IRS(I1,I2)=0  
    3   CONTINUE
    2 CONTINUE
      DO 4 I1=1,20
	IW1(1,I1)=0
	IW1(2,I1)=0
        DO 5 I2=1,20
          IRL(I1,I2)=0  
    5   CONTINUE
    4 CONTINUE
      DO 1 I=IP,IG
        KL=I-1
        IF(IA.EQ.IB) THEN
*
          IF(IOCASE.EQ.1) THEN
* Spin-spin operator
            CALL SSC(IG,KL,IC,ID,IA,IA,IA,IB,IC,ID,KL,TWO41)
          ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
            CALL SOOC(IG,KL,IC,ID,IA,IA,IA,IB,IC,ID,KL,TWO41)
          ENDIF
*
        ELSEIF(IC.EQ.ID) THEN
*
          IF(IOCASE.EQ.1) THEN
* Spin-spin operator
            CALL SSC(IG,KL,IA,IB,IC,IC,IA,IB,IC,ID,KL,TWO42)
          ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
            CALL SOOC(IG,KL,IA,IB,IC,IC,IA,IB,IC,ID,KL,TWO42)
          ENDIF
*
        ELSE
          WRITE(6,'(A)')  ' ERROR IN SUBROUTINE TWOPARTICLE4  '
          STOP
        ENDIF
    1 CONTINUE
      RETURN
      END
