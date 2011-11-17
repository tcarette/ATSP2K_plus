*
*     -------------------------------------------------------------
*      T W O P A R T I C L E 3 2
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
*                                              N'2 = N2 + 1        *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE TWOPARTICLE32(IA,IB,IIA,IIB,IIC,IID)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      COMMON/PERMAT/IRS(2,2),IRL(20,20),RS(2,2),RL(20,20)
      COMMON/KAMPAS/ IW1(2,20),IW2(2,20),IWAA(2,2,20,20),
     :RW1(2,20),RW2(2,20),RWAA(2,2,20,20)
      COMMON /CASEOP/ IOCASE
      EXTERNAL TWO32
      IF(IHSH.LE.1)RETURN
      CALL HIBFF(IA,IB,IA,IA,2)
      IF(ITTK(1,IK1(6),ID1(6)).EQ.0) RETURN
      IF(ITTK(2*IK1(3),IK1(5),ID1(5)).EQ.0) RETURN
      IF(IABS(IK2(6)-ID2(6)).GT.3) RETURN
      IF(IABS(IK2(5)-ID2(5)).GT.6*IK2(3)) RETURN
      IF(IOCASE.EQ.1) THEN
* Spin-Spin operator
        IP=ITREXG2(LJ(IA),LJ(IB),LJ(IB),LJ(IB),IKK)+1
      ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
        IP=ITREXG(LJ(IA),LJ(IB),LJ(IB),LJ(IB),IKK)+1
      ENDIF
      IF(IKK.LE.0)RETURN
      IG=IP+IKK-1
      DO 2 I1=1,2
        IRS(I1,1)=0
        DO 3 I2=1,20
          IRL(I2,1)=0
          DO 4 I3=1,2
            DO 5 I4=1,20
	      IWAA(I1,I3,I2,I4)=0
    5       CONTINUE
    4     CONTINUE
    3   CONTINUE
    2 CONTINUE
      DO 1 I=IP,IG
        KL=I-1
*
        IF(IOCASE.EQ.1) THEN
* Spin-spin operator
          CALL SSC(IG,KL,IA,IB,IB,IB,IB,IB,IB,IA,KL,TWO32)
        ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
          CALL SOOC(IG,KL,IA,IB,IB,IB,IB,IB,IB,IA,KL,TWO32)
        ENDIF
*
    1 CONTINUE
      RETURN
      END
