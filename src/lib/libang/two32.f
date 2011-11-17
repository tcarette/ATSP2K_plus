*
*     -------------------------------------------------------------
*      T W O 3 2
*     -------------------------------------------------------------
*                                                                  *
*                                                                  *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
*                                              N'2 = N2 + 1        *
*                                                                  *
*                                                                  *
*     CASES 2221   + + - -        TRANSFORM TO  1222   - + + -     *
*           2212                                1222               *
*                                                                  *
*                                                                  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE TWO32(KL1,KS1,KL2,KS2,K,IA,IB,INE1,INE2,INE3,CA1,CB1)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/PERMAT/IRS(2,2),IRL(20,20),RS(2,2),RL(20,20)
      COMMON /CASEOP/ IOCASE
      CA1=ZERO
      CB1=ZERO
      IAA=MIN0(IA,IB)
      IBB=MAX0(IA,IB)
*
      IF(IOCASE.EQ.1) THEN
* Spin-spin operator
        IPS1=2
      ELSE
        IPS1=1
      ENDIF
      DO 2 K2KS3=1,3,2
        ISES1=IXJTIK(2*KS1,2*KS2,2*K,1,K2KS3,1)
        ISES2=IXJTIK(2*KS2,2*KS1,2*K,1,K2KS3,1)
        IF(ISES1+ISES2.EQ.0) GO TO 2
        IF(IA.EQ.IAA) THEN
          KIS1=1
          KIS2=K2KS3
        ELSE
          KIS1=K2KS3
          KIS2=1
        ENDIF
        KIIS1=1
	IF(K2KS3.EQ.3)KIIS1=2
        IF(IRS(KIIS1,1).EQ.0) THEN
	  IRS(KIIS1,1)=1
          CALL RLSP2(3,IAA,IBB,KIS1,KIS2,2*K,0,IAT,REC)
          IF(IAT.EQ.0) THEN
            RECS=ZERO
          ELSE
            CALL RLSP2(3,IAA,IBB,KIS1,KIS2,2*K,1,IAT,RECS)
          ENDIF
	  RS(KIIS1,1)=RECS
        ELSE
          RECS=RS(KIIS1,1)
        ENDIF
        IF(DABS(RECS).LT.EPS) GO TO 2
        RECS=RECS*DSQRT(DBLE(K2KS3+1))
	IF(ISES1.NE.0) THEN
          CALL SIXJ(2*KS1,2*KS2,2*K,1,K2KS3,1,0,SN1)
        ELSE
	  SN1=ZERO
        ENDIF
	IF(ISES2.NE.0) THEN
          CALL SIXJ(2*KS2,2*KS1,2*K,1,K2KS3,1,0,SN2)
        ELSE
	  SN2=ZERO
        ENDIF
      DO 4 I1=IPS1,2
        KKS1=I1-1
        ISES3=IXJTIK(2*KKS1,1,1,2*KS1,K2KS3,1)
        ISES4=IXJTIK(2*KKS1,1,1,2*KS2,K2KS3,1)
        IF(ISES3+ISES4.EQ.0) GO TO 4
	IF(ISES1*ISES3.NE.0) THEN
          CALL SIXJ(2*KKS1,1,1,2*KS1,K2KS3,1,0,SN3)
          SN3=SN3*DSQRT(DBLE(2*KKS1+1))
          IF(MOD(KKS1+K2KS3,2).NE.0) SN3=-SN3
        ELSE
	  SN3=ZERO
        ENDIF
	IF(ISES2*ISES4.NE.0) THEN
          CALL SIXJ(2*KKS1,1,1,2*KS2,K2KS3,1,0,SN4)
          SN4=SN4*DSQRT(DBLE(2*KKS1+1))
          IF(MOD(2*KKS1+KS1+KS2+K2KS3,2).NE.0) SN4=-SN4
        ELSE
	  SN4=ZERO
        ENDIF
        IP1=ITREXG(K,LJ(IA),IK2(5)/2,ID2(5)/2,IKK)+1
        IF(IKK.LE.0)GO TO 4
        IG1=IP1+IKK-1
      DO 5 I5=IP1,IG1
        KKL3=I5-1
        ISES1L=IXJTIK(2*KL1,2*KL2,2*K,2*ID1(3),2*KKL3,2*ID2(3))
        ISES2L=IXJTIK(2*KL2,2*KL1,2*K,2*ID1(3),2*KKL3,2*ID2(3))
        IF(ISES1L+ISES2L.EQ.0) GO TO 5
        IF(IA.EQ.IAA) THEN
          KIL1=2*ID1(3)
          KIL2=2*KKL3
        ELSE
          KIL1=2*KKL3
          KIL2=2*ID1(3)
        ENDIF
        IF(IRL(I5,1).EQ.0) THEN
        CALL RLSP2(2,IAA,IBB,KIL1,KIL2,2*K,0,IAT,REC)
          IF(IAT.EQ.0) THEN
            IRL(I5,1)=1
            RECL=ZERO
	    RL(I5,1)=RECL
          ELSE
            RECL=ONE
          ENDIF
        ELSE
	  RECL=RL(I5,1)
        ENDIF
	IF(DABS(RECL).LT.EPS) GO TO 5
        IF(IRL(I5,1).EQ.0) THEN
          IRL(I5,1)=1
          CALL RLSP2(2,IAA,IBB,KIL1,KIL2,2*K,1,IAT,RECL)
	  RL(I5,1)=RECL
        ENDIF
	IF(DABS(RECL).LT.EPS) GO TO 5
        IP2=ITREXG(LJ(IB),LJ(IB),LJ(IB),KKL3,IKK)+1
        IF(IKK.LE.0) GO TO 5
        IG2=IP2+IKK-1
      DO 7 I3=IP2,IG2
        KKL1=I3-1
	IAT=1
        IF(IOCASE.EQ.1) THEN
* Spin-spin operator
          IF(MOD(KKL1+1,2).NE.0) IAT=0
        ENDIF
        IF(IAT.EQ.0) GO TO 7
        ISES3L=IXJTIK(2*KKL1,2*ID2(3),2*ID2(3),2*KL1,2*KKL3,2*ID2(3))
        ISES4L=IXJTIK(2*KKL1,2*ID2(3),2*ID2(3),2*KL2,2*KKL3,2*ID2(3))
        IF(ISES3L+ISES4L.EQ.0) GO TO 7
        BKKS2=HALF*DBLE(K2KS3)
        CALL WA1A2LSP(IAA,IBB,KKL1,KKS1,KKL3,BKKS2,
     :                                       HALF,HALF,-HALF,-HALF,W)
        IF(DABS(W).LT.EPS) GO TO 7
	IF(ISES1*ISES3*ISES1L*ISES3L.NE.0) THEN
          CALL SIXJ(2*KL1,2*KL2,2*K,2*ID1(3),2*KKL3,2*ID2(3),0,SLA1)
          CALL 
     :        SIXJ(2*KKL1,2*ID2(3),2*ID2(3),2*KL1,2*KKL3,2*ID2(3),0,S3)
          SLA1=SLA1*S3*DSQRT(DBLE(2*KKL1+1))
          IF(MOD(KKL1+ID1(3)+ID2(3)+2*KKL3,2).NE.0) SLA1=-SLA1
        ELSE
	  SLA1=ZERO
        ENDIF
	IF(ISES2*ISES4*ISES2L*ISES4L.NE.0) THEN
          CALL SIXJ(2*KL2,2*KL1,2*K,2*ID1(3),2*KKL3,2*ID2(3),0,SLB1)
          CALL 
     :        SIXJ(2*KKL1,2*ID2(3),2*ID2(3),2*KL2,2*KKL3,2*ID2(3),0,S4)
          SLB1=SLB1*S4*DSQRT(DBLE(2*KKL1+1))
          IF(MOD(2*KKL1+KL1+KL2+2*KKL3+ID1(3)+ID2(3),2).NE.0)SLB1=-SLB1
        ELSE
	  SLB1=ZERO
        ENDIF
        CA=DSQRT(DBLE(2*KKL3+1))*SN1*SN3*RECS*SLA1*W*RECL
        CB=DSQRT(DBLE(2*KKL3+1))*SN2*SN4*RECS*SLB1*W*RECL
        IF(IA.EQ.IAA) THEN
          IF(MOD(KIL1+KIS1+KIL2+KIS2-2*K-2*K+2,4).NE.0) CA=-CA
          IF(MOD(KIL1+KIS1+KIL2+KIS2-2*K-2*K+2,4).NE.0) CB=-CB
        ENDIF
        CA1=CA1+CA
        CB1=CB1+CB
    7 CONTINUE
    5 CONTINUE
    4 CONTINUE
    2 CONTINUE
      CA1=CA1*HALF
      CB1=CB1*HALF
      RETURN
      END
