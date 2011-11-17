*
*     -------------------------------------------------------------
*      T W O 5 5 
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
*                                              N'2 = N2 + 1        *
*     ( IREZ = 1 )                             N'3 = N3 - 1        *
*                                              N'4 = N4 + 1        *
*                                                                  *
*     CASES 1423   + + - -        TRANSFORM TO  1234     + - - +   *
*           4132                                1234               *
*                                                                  *
*                                                                  *
*                                                                  *
*     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 + 1        *
*                                              N'2 = N2 - 1        *
*     ( IREZ = 2 )                             N'3 = N3 + 1        *
*                                              N'4 = N4 - 1        *
*                                                                  *
*     CASES 2314   + + - -        TRANSFORM TO  1234     - + + -   *
*           3241                                1234               *
*                                                                  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE TWO55(KL1,KS1,KL2,KS2,K,IA,IB,IC,ID,IREZ,CA1,CB1)
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
        IF(IW1(KS1+1,3).EQ.0) THEN
	  IW1(KS1+1,3)=1
          CALL RLSP4A(3,IA,IB,IC,1,1,2*KS1,0,IAT,R)
          IF(IAT.EQ.0) THEN
	    RECS2A=ZERO
          ELSE
            CALL RLSP4A(3,IA,IB,IC,1,1,2*KS1,1,IAT,RECS2A)
	  ENDIF
	  RW1(KS1+1,3)=RECS2A
        ELSE
          RECS2A=RW1(KS1+1,3)
        ENDIF
        IF(IW1(KS2+1,3).EQ.0) THEN
	  IW1(KS2+1,3)=1
          CALL RLSP4A(3,IA,IB,IC,1,1,2*KS2,0,IAT,R)
          IF(IAT.EQ.0) THEN
	    RECS2B=ZERO
          ELSE
            CALL RLSP4A(3,IA,IB,IC,1,1,2*KS2,1,IAT,RECS2B)
	  ENDIF
	  RW1(KS2+1,3)=RECS2B
        ELSE
          RECS2B=RW1(KS2+1,3)
        ENDIF
        ISES3=1
        ISES4=1
        IF(DABS(RECS2A).LT.EPS) ISES3=0
        IF(DABS(RECS2B).LT.EPS) ISES4=0
        IF(ISES3+ISES4.EQ.0) GO TO 3
        ISES1=0
        ISES2=0
        DO 1 K2KS3=IPS1,IGS1,2
          IF(ISES3.NE.0) ISES1=IXJTIK(2*KS2,2*KS1,2*K,K2KS3,1,1)
          IF(ISES4.NE.0) ISES2=IXJTIK(2*KS1,2*KS2,2*K,K2KS3,1,1)
          IF(ISES1+ISES2.NE.0) THEN
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
            IF(DABS(RECS1).GT.EPS) THEN
              IF(ISES1.NE.0) THEN
                IF(IRS(KS1+1,KIIS1).EQ.0) THEN
                  IRS(KS1+1,KIIS1)=1
                  IAT=0
                  CALL DLSA4(3,IB,IC,2*KS1,1,K2KS3,0,IAT,R)
                  IF(IAT.EQ.0) THEN
	            RECS3A=ZERO
                  ELSE
                    CALL DLSA4(3,IB,IC,2*KS1,1,K2KS3,1,IAT,RECS3A)
                  ENDIF
	          RS(KS1+1,KIIS1)=RECS3A
                ELSE
                  RECS3A=RS(KS1+1,KIIS1)
                ENDIF
	        IF(DABS(RECS3A).LT.EPS) ISES1=0
              ENDIF
              IF(ISES2.NE.0) THEN
                IF(IRS(KS2+1,KIIS1).EQ.0) THEN
                  IRS(KS2+1,KIIS1)=1
                  IAT=0
                  CALL DLSA4(3,IB,IC,2*KS2,1,K2KS3,0,IAT,R)
                  IF(IAT.EQ.0) THEN
	            RECS3B=ZERO
                  ELSE
                    CALL DLSA4(3,IB,IC,2*KS2,1,K2KS3,1,IAT,RECS3B)
                  ENDIF
	          RS(KS2+1,KIIS1)=RECS3B
                ELSE
                  RECS3B=RS(KS2+1,KIIS1)
                ENDIF
	        IF(DABS(RECS3B).LT.EPS) ISES2=0
              ENDIF
            ENDIF
            IF(ISES1+ISES2.NE.0) THEN
              IF(ISES1.NE.0) THEN
                CALL SIXJ(2*KS2,2*KS1,2*K,K2KS3,1,1,0,SN1)
                SN1=SN1*DSQRT(DBLE(K2KS3+1))
                CAS1=CAS1+SN1*RECS1*RECS2A*RECS3A
              ENDIF
              IF(ISES2.NE.0) THEN
                CALL SIXJ(2*KS1,2*KS2,2*K,K2KS3,1,1,0,SN2)
                SN2=SN2*DSQRT(DBLE(K2KS3+1))
                CBS1=CBS1+SN2*RECS1*RECS2B*RECS3B
              ENDIF
            ENDIF
          ENDIF
    1   CONTINUE
    3   CONTINUE
	CAS1=CAS1/DSQRT(DBLE(2*KS1+1))
	CBS1=CBS1/DSQRT(DBLE(2*KS2+1))
        IF(IREZ.EQ.1) THEN
          IF(MOD(KS1+KS2,2).NE.0) CAS1=-CAS1
        ELSE
          IF(MOD(KS1+KS2,2).NE.0) CBS1=-CBS1
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
        IF(IW2(2,KL1+1).EQ.0) THEN
	  IW2(2,KL1+1)=1
          CALL RLSP4A(2,IA,IB,IC,2*ID1(3),2*ID2(3),2*KL1,0,IAT,R)
          IF(IAT.EQ.0) THEN
	    RECL2A=ZERO
          ELSE
            CALL RLSP4A
     :               (2,IA,IB,IC,2*ID1(3),2*ID2(3),2*KL1,1,IAT,RECL2A)
	  ENDIF
	  RW2(2,KL1+1)=RECL2A
        ELSE
          RECL2A=RW2(2,KL1+1)
        ENDIF
        IF(IW2(2,KL2+1).EQ.0) THEN
	  IW2(2,KL2+1)=1
          CALL RLSP4A(2,IA,IB,IC,2*ID1(3),2*ID2(3),2*KL2,0,IAT,R)
          IF(IAT.EQ.0) THEN
	    RECL2B=ZERO
          ELSE
            CALL RLSP4A
     :               (2,IA,IB,IC,2*ID1(3),2*ID2(3),2*KL2,1,IAT,RECL2B)
	  ENDIF
	  RW2(2,KL2+1)=RECL2B
        ELSE
          RECL2B=RW2(2,KL2+1)
        ENDIF
        ISES3=1
        ISES4=1
        IF(DABS(RECL2A).LT.EPS) ISES3=0
        IF(DABS(RECL2B).LT.EPS) ISES4=0
        IF(ISES3+ISES4.EQ.0) GO TO 4
        ISES1L=0
        ISES2L=0
        DO 2 I3=IP1,IG1
          KKL3=I3-1
          IF(ISES3.NE.0) 
     :         ISES1L=IXJTIK(2*KL2,2*KL1,2*K,2*KKL3,2*ID4(3),2*ID3(3))
          IF(ISES4.NE.0) 
     :         ISES2L=IXJTIK(2*KL1,2*KL2,2*K,2*KKL3,2*ID4(3),2*ID3(3))
          IF(ISES1L+ISES2L.NE.0) THEN
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
            IF(DABS(RECL1).GT.EPS) THEN
              IF(ISES1L.NE.0) THEN
                IF(IRL(KL1+1,I3).EQ.0) THEN
                  IRL(KL1+1,I3)=1
                  IAT=0
                  CALL DLSA4(2,IB,IC,2*KL1,2*ID3(3),2*KKL3,0,IAT,R)
                  IF(IAT.EQ.0) THEN
	            RECL3A=ZERO
                  ELSE
                    CALL DLSA4
     :                   (2,IB,IC,2*KL1,2*ID3(3),2*KKL3,1,IAT,RECL3A)
                  ENDIF
	          RL(KL1+1,I3)=RECL3A
                ELSE
                  RECL3A=RL(KL1+1,I3)
                ENDIF
	        IF(DABS(RECL3A).LT.EPS) ISES1L=0
              ELSE
                ISES1L=0
              ENDIF
              IF(ISES2L.NE.0) THEN
                IF(IRL(KL2+1,I3).EQ.0) THEN
                  IRL(KL2+1,I3)=1
                  IAT=0
                  CALL DLSA4(2,IB,IC,2*KL2,2*ID3(3),2*KKL3,0,IAT,R)
                  IF(IAT.EQ.0) THEN
	            RECL3B=ZERO
                  ELSE
                    CALL DLSA4
     :                   (2,IB,IC,2*KL2,2*ID3(3),2*KKL3,1,IAT,RECL3B)
                  ENDIF
	          RL(KL2+1,I3)=RECL3B
                ELSE
                  RECL3B=RL(KL2+1,I3)
                ENDIF
	        IF(DABS(RECL3B).LT.EPS) ISES2L=0
              ELSE
                ISES2L=0
              ENDIF
              IF(ISES1L+ISES2L.NE.0) THEN
	        IF(ISES1L.NE.0) THEN
                  CALL SIXJ
     :                  (2*KL2,2*KL1,2*K,2*KKL3,2*ID4(3),2*ID3(3),0,S3)
                  CA1=CA1+S3*RECL1*RECL2A*RECL3A*DSQRT(DBLE(2*KKL3+1))
                ENDIF
	        IF(ISES2L.NE.0) THEN
                  CALL SIXJ
     :                  (2*KL1,2*KL2,2*K,2*KKL3,2*ID4(3),2*ID3(3),0,S4)
                  CB1=CB1+S4*RECL1*RECL2B*RECL3B*DSQRT(DBLE(2*KKL3+1))
                ENDIF
              ENDIF
            ENDIF
          ENDIF
    2   CONTINUE
    4   CONTINUE
	CA1=CA1/DSQRT(DBLE(2*KL1+1))
	CB1=CB1/DSQRT(DBLE(2*KL2+1))
        IF(IREZ.EQ.1) THEN
          IF(MOD(KL1+KL2+1,2).NE.0) CA1=-CA1
          CB1=-CB1
        ELSE
          IF(MOD(ID1(3)+ID2(3)+ID3(3)+ID4(3)+1,2).NE.0) CA1=-CA1
          IF(MOD(ID1(3)+ID2(3)+ID3(3)+ID4(3)+KL1+KL2+1,2).NE.0)CB1=-CB1
        ENDIF
        RWAA(1,1,KL1+1,KL2+1)=CA1
        RWAA(1,2,KL1+1,KL2+1)=CB1
      ELSE
        CA1=RWAA(1,1,KL1+1,KL2+1)
        CB1=RWAA(1,2,KL1+1,KL2+1)
      ENDIF
      CA1=HALF*CA1*CAS1*RW1(1,1)
      CB1=HALF*CB1*CBS1*RW1(1,1)
      RETURN
      END
