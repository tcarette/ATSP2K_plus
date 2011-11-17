*
*     -------------------------------------------------------------
*      T W O 2
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF ONE PARTICLE OPERATOR IN CASE :       N'1 = N1 - 2        *
*                                              N'2 = N2 + 2        *
*                                                                  *
*                                                     1122         *
*                                                                  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE TWO2(KL1,KS1,KL2,KS2,K,IA,IB,INE1,INE2,INE3,C1,C2)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/PERMAT/IRS(2,2),IRL(20,20),RS(2,2),RL(20,20)
      COMMON /CASEOP/ IOCASE
      C1=ZERO
      C2=ZERO
      IF(IHSH.LE.1) RETURN
      IF(IA.EQ.IB) RETURN
      IP1=ITREXG(LJ(IA),LJ(IA),IK1(5)/2,ID1(5)/2,IKK)+1
      IF(IKK.LE.0) RETURN
      IG1=IP1+IKK-1
      IP2=ITREXG(LJ(IB),LJ(IB),IK2(5)/2,ID2(5)/2,IKK)+1
      IF(IKK.LE.0) RETURN
      IG2=IP2+IKK-1
      IIA=MIN0(IA,IB)
      IIB=MAX0(IA,IB)
      IF(IOCASE.EQ.1) THEN
* Spin - spin      operator
        IPS1=2
        IPS2=2
      ELSE
* Spin-other-orbit operator
        IPS1=1
        IPS2=1
      ENDIF
      DO 2 I1=IPS1,2
      KKS1=I1-1
      DO 3 I2=IPS2,2
        KKS2=I2-1
        IF(IA.NE.IIA) THEN
          KIS1=KKS2
          KIS2=KKS1
        ELSE
          KIS1=KKS1
          KIS2=KKS2
        ENDIF
        CALL NINESS(KS1,KS2,KKS1,KKS2,K,SN1)
        IF(DABS(SN1).LT.EPS) GO TO 3
        IF(IRS(KIS1+1,KIS2+1).EQ.0) THEN
	  IRS(KIS1+1,KIS2+1)=1
          CALL RLSP2(3,IIA,IIB,2*KIS1,2*KIS2,2*K,0,IAT,REC)
          IF(IAT.EQ.0) THEN
	    RS(KIS1+1,KIS2+1)=ZERO
            RECS=ZERO
          ELSE
            CALL RLSP2(3,IIA,IIB,2*KIS1,2*KIS2,2*K,1,IAT,RECS)
	    RS(KIS1+1,KIS2+1)=RECS
          ENDIF
        ELSE
	  RECS=RS(KIS1+1,KIS2+1)
        ENDIF
        IF(DABS(RECS).LT.EPS) GO TO 3
      DO 4 I3=IP1,IG1
      KKL1=I3-1
      DO 5 I4=IP2,IG2
        KKL2=I4-1
	IAT=1
        IF(IOCASE.EQ.1) THEN
* Spin-spin operator
          IF(MOD(KKL1+KKL2,2).NE.0) IAT=0
        ENDIF
        IF(IAT.EQ.0) GO TO 5
        IF(KL1.EQ.KL2) THEN
          IF(MOD(KKL1+KKL2+K,2).NE.0) IAT=0
        ENDIF
        IF(IAT.EQ.0) GO TO 5
        IF(IA.NE.IIA) THEN
          KIL1=KKL2
          KIL2=KKL1
        ELSE
          KIL1=KKL1
          KIL2=KKL2
        ENDIF
        CALL NINELS(2*ID1(3),2*ID2(3),2*KL1,2*ID1(3),2*ID2(3),
     :                              2*KL2,2*KKL1,2*KKL2,2*K,1,IN,SN2)
        IF(IN.EQ.0) GO TO 5
        IF(IRL(KIL1+1,KIL2+1).EQ.0) THEN
	  IRL(KIL1+1,KIL2+1)=1
          CALL RLSP2(2,IIA,IIB,2*KIL1,2*KIL2,2*K,0,IAT,REC)
          IF(IAT.EQ.0) THEN
            RECL=ZERO
          ELSE
            CALL RLSP2(2,IIA,IIB,2*KIL1,2*KIL2,2*K,1,IAT,RECL)
          ENDIF
	  RL(KIL1+1,KIL2+1)=RECL
        ELSE
	  RECL=RL(KIL1+1,KIL2+1)
        ENDIF
        IF(DABS(RECL).LT.EPS) GO TO 5
        CALL W1W2LSP(KKL1,KKS1,KKL2,KKS2,HALF,HALF,-HALF,-HALF,W)
        IF(DABS(W).LT.EPS) GO TO 5
        W=W*DSQRT(DBLE((2*KKL1+1)*(2*KKL2+1)*(2*KKS1+1)*(2*KKS2+1)))
        CALL NINELS(2*ID1(3),2*ID2(3),2*KL1,2*ID1(3),2*ID2(3),
     :                            2*KL2,2*KKL1,2*KKL2,2*K,0,IN,SN2)
        C=-RECS*SN1*W*RECL*SN2*HALF
        IF(IA.NE.IIA) THEN
          IF(MOD(KKL1+KKS1+KKL2+KKS2-K-K,2).NE.0) C=-C
        ENDIF
        C1=C1+C
    5 CONTINUE
    4 CONTINUE
    3 CONTINUE
    2 CONTINUE
      RETURN
      END
