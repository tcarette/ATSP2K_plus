*
*     -------------------------------------------------------------
*      T W O 4 1
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 + 1        *
*                                              N'2 = N2 + 1        *
*                                              N'3 = N3 - 2,       *
*                                                                  *
*     CASES 3312   + + - -        TRANSFORM TO  1233   - - + +     *
*           3321                                1233               *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE TWO41(KL1,KS1,KL2,KS2,K,IA,IB,IC,INE1,INE2,CA1,CB1)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      COMMON/TRK2/BD3(3),BD4(3),BK3(3),BK4(3),
     *ID3(7),ID4(7),IK3(7),IK4(7)
      COMMON /CASEOP/ IOCASE
      COMMON/PERMAT/IRS(2,2),IRL(20,20),RS(2,2),RL(20,20)
      CA1=ZERO
      CB1=ZERO
      IF(IHSH.LE.2)RETURN
      IP1=ITREXG(LJ(IC),LJ(IC),IK3(5)/2,ID3(5)/2,IKK)+1
      IF(IKK.LE.0) RETURN
      IG1=IP1+IKK-1
*
      IF(IOCASE.EQ.1) THEN
* Spin-spin operator
        IPS1=2
        IPS2=2
      ELSE
        IPS1=1
        IPS2=1
      ENDIF
      DO 2 I1=IPS1,2
      KKS1=I1-1
      DO 3 I2=IPS2,2
        KKS2=I2-1
        CALL NINESS(KS1,KS2,KKS2,KKS1,K,SN1)
        IF(DABS(SN1).LT.EPS) GO TO 3
        SN1=SN1*DSQRT(DBLE((2*KKS1+1)*(2*KKS2+1)))
        IF(IRS(I1,I2).EQ.0) THEN
	  IRS(I1,I2)=1
          CALL RLSP3(3,IA,IB,IC,1,1,2*KKS1,2*KKS2,2*K,RECS)
	  RS(I1,I2)=RECS
        ELSE
	  RECS=RS(I1,I2)
        ENDIF
        IF(DABS(RECS).LT.EPS) GO TO 3
        IP=IABS(LJ(IB)-LJ(IA))+1
        IG=LJ(IB)+LJ(IA)+1
      DO 4 I3=IP1,IG1
        KKL2=I3-1
        CALL A1A2W3LSP(IC,IA,IB,IC,KKL2,KKS2,-HALF,-HALF,HALF,HALF,W)
        IF(DABS(W).LT.EPS) GO TO 4 
      DO 5 I4=IP,IG
        KKL1=I4-1
	IAT=1
        IF(IOCASE.EQ.1) THEN
* Spin-spin operator
          IF(MOD(1+KKL2,2).NE.0) IAT=0
        ENDIF
        IF(IAT.EQ.0) GO TO 5
        CALL NINELS(2*ID3(3),2*ID1(3),2*KL1,2*ID3(3),2*ID2(3),2*KL2,
     :         2*KKL2,2*KKL1,2*K,1,INA,S)
        CALL NINELS(2*ID3(3),2*ID2(3),2*KL1,2*ID3(3),2*ID1(3),2*KL2,
     :         2*KKL2,2*KKL1,2*K,1,INB,S)
        IF((INA+INB).EQ.0) GO TO 5
        IF(IRL(I4,I3).EQ.0) THEN
	  IRL(I4,I3)=1
          CALL 
     :       RLSP3(2,IA,IB,IC,2*ID1(3),2*ID2(3),2*KKL1,2*KKL2,2*K,RECL)
	  RL(I4,I3)=RECL
        ELSE
	  RECL=RL(I4,I3)
        ENDIF
        IF(DABS(RECL).LT.EPS) GO TO 5
        IF(INA.NE.0) THEN
        CALL NINELS(2*ID3(3),2*ID1(3),2*KL1,2*ID3(3),2*ID2(3),2*KL2,
     :         2*KKL2,2*KKL1,2*K,0,INA,SLA1)
          SLA1=SLA1*DSQRT(DBLE((2*KKL1+1)*(2*KKL2+1)))
          IF(MOD(KKL1+KKS1+KKL2+KKS2,2).NE.0) SLA1=-SLA1
        ELSE
          SLA1=ZERO
        ENDIF
        IF(INB.NE.0) THEN
        CALL NINELS(2*ID3(3),2*ID2(3),2*KL1,2*ID3(3),2*ID1(3),2*KL2,
     :         2*KKL2,2*KKL1,2*K,0,INB,SLB1)
          SLB1=SLB1*DSQRT(DBLE((2*KKL1+1)*(2*KKL2+1)))
          IF(MOD(ID1(3)+ID2(3)+KKL2+KKS2,2).NE.0) SLB1=-SLB1
        ELSE
          SLB1=ZERO
        ENDIF
        CA1=CA1+SN1*RECS*SLA1*W*RECL
        CB1=CB1+SN1*RECS*SLB1*W*RECL
    5 CONTINUE
    4 CONTINUE
    3 CONTINUE
    2 CONTINUE
      CA1=-CA1*HALF
      CB1=-CB1*HALF
      RETURN
      END
