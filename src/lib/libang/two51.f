*
*     -------------------------------------------------------------
*      T W O 5 1 
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
*                                              N'2 = N2 - 1        *
*     ( IREZ = 1 )                             N'3 = N3 + 1        *
*                                              N'4 = N4 + 1        *
*                                                                  *
*     CASES 1234   + + - -        TRANSFORM TO  1234     + + - -   *
*           2143                                1234               *
*                                                                  *
*                                                                  *
*                                                                  *
*     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 + 1        *
*                                              N'2 = N2 + 1        *
*     ( IREZ = 2 )                             N'3 = N3 - 1        *
*                                              N'4 = N4 - 1        *
*                                                                  *
*     CASES 3412   + + - -        TRANSFORM TO  1234     - - + +   *
*           4321                                1234               *
*                                                                  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE TWO51(KL1,KS1,KL2,KS2,K,IA,IB,IC,ID,IREZ,CA1,CB1)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      COMMON/TRK2/BD3(3),BD4(3),BK3(3),BK4(3),
     *ID3(7),ID4(7),IK3(7),IK4(7)
      COMMON/KAMPAS/ IW1(2,20),IW2(2,20),IWAA(2,2,20,20),
     :RW1(2,20),RW2(2,20),RWAA(2,2,20,20)
      COMMON/PERMAT/IRS(2,2),IRL(20,20),RS(2,2),RL(20,20)
      COMMON /CASEOP/ IOCASE
      CA1=ZERO
      CB1=ZERO
      IF(IHSH.LE.3)RETURN
      JN=IHSH+IC-1
      JI=J1QN1(JN,3)-1
      JIS=J1QN2(JN,3)-1
      IPS1=ITREXG(1,2*K,JI,JIS,IKK1)
      IF(IKK1.EQ.0) RETURN
      IF(IPS1.GT.3) RETURN
      IGS1=IPS1+IKK1-1
      IF(IGS1.GT.3) IGS1=3
      IF(IOCASE.EQ.1) THEN
* Spin-spin operator
        IF(IGS1.NE.3) RETURN
        IPS2=2
      ELSE
        IPS2=1
      ENDIF
      JI=(J1QN1(JN,2)-1)/2
      JIS=(J1QN2(JN,2)-1)/2
      IP1=ITREXG(LJ(ID),K,JI,JIS,IKK1)+1
      IF(IKK1.EQ.0) RETURN
      IG1=IP1+IKK1-1
      IF(IW1(KS1+1,KS2+4).EQ.0) THEN
        IW1(KS1+1,KS2+4)=1
        CAS1=ZERO
        CBS1=ZERO
        DO 1 K2KS3=IPS1,IGS1,2
          ISES1=IXJTIK(2*KS2,2*KS1,2*K,K2KS3,1,1)
          ISES2=IXJTIK(2*KS1,2*KS2,2*K,K2KS3,1,1)
          IF((ISES1+ISES2).EQ.0) GO TO 1
          KIIS1=1
          IF(K2KS3.EQ.3)KIIS1=2
          IF(IW1(KIIS1,2).EQ.0) THEN
	    IW1(KIIS1,2)=1
            CALL RLSP4B(3,IC,ID,K2KS3,1,2*K,0,IAT,R)
            IF(IAT.EQ.0) THEN
	      RECS1=ZERO
            ELSE
              CALL RLSP4B(3,IC,ID,K2KS3,1,2*K,1,IAT,RECS1)
	    ENDIF
	    RW1(KIIS1,2)=RECS1
          ELSE
            RECS1=RW1(KIIS1,2)
          ENDIF
          IF(DABS(RECS1).LT.EPS) GO TO 1
          IF(ISES1.NE.0) CALL SIXJ(2*KS2,2*KS1,2*K,K2KS3,1,1,0,SN1)
          IF(ISES2.NE.0) CALL SIXJ(2*KS1,2*KS2,2*K,K2KS3,1,1,0,SN2)
        DO 2 I1=IPS2,2
          KKS1=I1-1
          ISES3=IXJTIK(2*KKS1,1,1,2*KS1,1,K2KS3)
          ISES4=IXJTIK(2*KKS1,1,1,2*KS2,1,K2KS3)
          IF((ISES3+ISES4).NE.0) THEN
            IF(IRS(I1,KIIS1).EQ.0) THEN
              IRS(I1,KIIS1)=1
              IAT=0
              CALL DLSA4(3,IB,IC,2*KKS1,1,K2KS3,0,IAT,R)
              IF(IAT.EQ.0) THEN
	        RECS3=ZERO
              ELSE
                CALL DLSA4(3,IB,IC,2*KKS1,1,K2KS3,1,IAT,RECS3)
              ENDIF
	      RS(I1,KIIS1)=RECS3
            ELSE
              RECS3=RS(I1,KIIS1)
            ENDIF
            IF(DABS(RECS3).GT.EPS) THEN
              IF(IW1(I1,3).EQ.0) THEN
	        IW1(I1,3)=1
                CALL RLSP4A(3,IA,IB,IC,1,1,2*KKS1,0,IAT,R)
                IF(IAT.EQ.0) THEN
	          RECS2=ZERO
                ELSE
                  CALL RLSP4A(3,IA,IB,IC,1,1,2*KKS1,1,IAT,RECS2)
	        ENDIF
	        RW1(I1,3)=RECS2
              ELSE
                RECS2=RW1(I1,3)
              ENDIF
              IF(DABS(RECS2).GT.EPS) THEN
                IF(ISES1*ISES3.NE.0) THEN
                  CALL SIXJ(2*KKS1,1,1,2*KS1,1,K2KS3,0,SN3)
                  SN3=SN3*DSQRT(DBLE((K2KS3+1)*(2*KKS1+1)))
	          IF(MOD(KKS1,2).NE.0) SN3=-SN3
                  CAS1=CAS1+SN1*SN3*RECS1*RECS2*RECS3
                ENDIF
                IF(ISES2*ISES4.NE.0) THEN
                  CALL SIXJ(2*KKS1,1,1,2*KS2,1,K2KS3,0,SN4)
                  SN4=SN4*DSQRT(DBLE((K2KS3+1)*(2*KKS1+1)))
	          IF(MOD(KKS1,2).NE.0) SN4=-SN4
                  CBS1=CBS1+SN2*SN4*RECS1*RECS2*RECS3
                ENDIF
              ENDIF
            ENDIF
          ENDIF
    2   CONTINUE
    1   CONTINUE
	IF(IREZ.EQ.1) THEN
          IF(MOD(KS1+KS2,2).NE.0) CBS1=-CBS1
        ELSE
	  IF(MOD(KS1+KS2,2).NE.0) CAS1=-CAS1
        ENDIF
        RW1(KS1+1,KS2+4)=CAS1
        RW1(KS1+1,KS2+6)=CBS1
      ELSE
        CAS1=RW1(KS1+1,KS2+4)
        CBS1=RW1(KS1+1,KS2+6)
      ENDIF
      ISES1=1 
      ISES2=1 
      IF(DABS(CAS1).LT.EPS) ISES1=0 
      IF(DABS(CBS1).LT.EPS) ISES2=0 
      IF((ISES1+ISES2).EQ.0) RETURN
      IF(IWAA(1,1,KL1+1,KL2+1).EQ.0) THEN
        IWAA(1,1,KL1+1,KL2+1)=1
        DO 3 I3=IP1,IG1
          KKL3=I3-1
          IP2=ITREXG(LJ(IA),LJ(IB),KKL3,LJ(IC),IKK2)+1
          IF(IKK2.EQ.0) GO TO 3
          IG2=IP2+IKK2-1
          ISES1L=IXJTIK(2*KL2,2*KL1,2*K,2*KKL3,2*ID4(3),2*ID2(3))
          ISES2L=IXJTIK(2*KL1,2*KL2,2*K,2*KKL3,2*ID4(3),2*ID2(3))
          IF((ISES1L+ISES2L).EQ.0) GO TO 3
          IF(IW2(1,I3).EQ.0) THEN
	    IW2(1,I3)=1
            CALL RLSP4B(2,IC,ID,2*KKL3,2*ID4(3),2*K,0,IAT,R)
            IF(IAT.EQ.0) THEN
	      RECL1=ZERO
            ELSE
              CALL RLSP4B(2,IC,ID,2*KKL3,2*ID4(3),2*K,1,IAT,RECL1)
	    ENDIF
	    RW2(1,I3)=RECL1
          ELSE
            RECL1=RW2(1,I3)
          ENDIF
          IF(DABS(RECL1).LT.EPS) GO TO 3
	  IF(ISES1L.NE.0) 
     :       CALL SIXJ(2*KL2,2*KL1,2*K,2*KKL3,2*ID4(3),2*ID2(3),0,S3)
	  IF(ISES2L.NE.0) 
     :       CALL SIXJ(2*KL1,2*KL2,2*K,2*KKL3,2*ID4(3),2*ID2(3),0,S4)
        DO 4 I4=IP2,IG2
          KKL1=I4-1
          ISES3L=IXJTIK(2*KL1,2*ID3(3),2*ID1(3),2*KKL1,2*ID2(3),2*KKL3)
          ISES4L=IXJTIK(2*KL2,2*ID3(3),2*ID1(3),2*KKL1,2*ID2(3),2*KKL3)
          IF((ISES3L+ISES4L).NE.0) THEN
            IF(IRL(I4,I3).EQ.0) THEN
              IRL(I4,I3)=1
              IAT=0
              CALL DLSA4(2,IB,IC,2*KKL1,2*ID3(3),2*KKL3,0,IAT,R)
              IF(IAT.EQ.0) THEN
	        RECL3=ZERO
              ELSE
                CALL DLSA4(2,IB,IC,2*KKL1,2*ID3(3),2*KKL3,1,IAT,RECL3)
              ENDIF
	      RL(I4,I3)=RECL3
            ELSE
              RECL3=RL(I4,I3)
            ENDIF
            IF(DABS(RECL3).GT.EPS) THEN
              IF(IW2(2,I4).EQ.0) THEN
	        IW2(2,I4)=1
                CALL RLSP4A
     :                    (2,IA,IB,IC,2*ID1(3),2*ID2(3),2*KKL1,0,IAT,R)
                IF(IAT.EQ.0) THEN
	          RECL2=ZERO
                ELSE
                  CALL RLSP4A
     :              (2,IA,IB,IC,2*ID1(3),2*ID2(3),2*KKL1,1,IAT,RECL2)
	        ENDIF
	        RW2(2,I4)=RECL2
              ELSE
                RECL2=RW2(2,I4)
              ENDIF
              IF(DABS(RECL2).GT.EPS) THEN
	        IF(ISES1L*ISES3L.NE.0) THEN
                  CALL SIXJ
     :          (2*KL1,2*ID3(3),2*ID1(3),2*KKL1,2*ID2(3),2*KKL3,0,SLA1)
                  SLA1=SLA1*DSQRT(DBLE((2*KKL1+1)*(2*KKL3+1)))
	          IF(MOD(KKL1,2).NE.0) SLA1=-SLA1
                  CA1=CA1+S3*SLA1*RECL1*RECL2*RECL3
                ENDIF
	        IF(ISES2L*ISES4L.NE.0) THEN
                  CALL SIXJ
     :          (2*KL2,2*ID3(3),2*ID1(3),2*KKL1,2*ID2(3),2*KKL3,0,SLB1)
                  SLB1=SLB1*DSQRT(DBLE((2*KKL1+1)*(2*KKL3+1)))
                  IF(MOD(KKL1,2).NE.0) SLB1=-SLB1
                  CB1=CB1+S4*SLB1*RECL1*RECL2*RECL3
                ENDIF
              ENDIF
            ENDIF
          ENDIF
    4   CONTINUE
    3   CONTINUE
	IF(IREZ.EQ.1) THEN
	  IF(MOD(ID3(3)+ID4(3),2).NE.0) CA1=-CA1
          IF(MOD(ID3(3)+ID4(3)+KL1+KL2,2).NE.0) CB1=-CB1
        ELSE
          IF(MOD(ID1(3)+ID2(3)+KL1+KL2,2).NE.0) CA1=-CA1
	  IF(MOD(ID1(3)+ID2(3),2).NE.0) CB1=-CB1
        ENDIF
        RWAA(1,1,KL1+1,KL2+1)=CA1
        RWAA(1,2,KL1+1,KL2+1)=CB1
      ELSE
        CA1=RWAA(1,1,KL1+1,KL2+1)
        CB1=RWAA(1,2,KL1+1,KL2+1)
      ENDIF
      CA1=-CA1*CAS1*HALF*RW1(1,1)
      CB1=-CB1*CBS1*HALF*RW1(1,1)
      RETURN
      END
