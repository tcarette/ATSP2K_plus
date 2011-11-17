*                                                                  *
*     -------------------------------------------------------------
*      T W O 3 3 A
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
*                                              N'2 = N2 + 1        *
*                                                                  *
*                                                                  *
*     CASES 2313   + + - -        TRANSFORM TO  2133   + - + -     *
*           3231                                2133               *
*                                                                  *
*                                                                  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE TWO33A(KL1,KS1,KL2,KS2,K,IA,IB,IC,INE1,INE2,CA1,CB1)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      COMMON/TRK2/BD3(3),BD4(3),BK3(3),BK4(3),
     *ID3(7),ID4(7),IK3(7),IK4(7)
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/PERMAT/IRS(2,2),IRL(20,20),RS(2,2),RL(20,20)
      COMMON /CASEOP/ IOCASE
      CA1=ZERO
      CB1=ZERO
C
C     CASES 2313   + + - -        TRANSFORM TO  2133   + - + - 
C
      CALL A1A2W3LSP(IB,IA,IC,IC,KL2,KS2,-HALF,HALF,HALF,-HALF,W)
      IF(DABS(W).GT.EPS) THEN
        IF(IRS(KS1+1,KS2+1).EQ.0) THEN
          IRS(KS1+1,KS2+1)=1
          CALL RLSP3(3,IB,IA,IC,1,1,2*KS1,2*KS2,2*K,RES)
          RS(KS1+1,KS2+1)=RES
        ELSE
          RES=RS(KS1+1,KS2+1)
        ENDIF
        IF(DABS(RES).GT.EPS)  THEN
          IF(IRL(KL1+1,KL2+1).EQ.0) THEN
            IRL(KL1+1,KL2+1)=1
            CALL 
     :         RLSP3(2,IB,IA,IC,2*ID2(3),2*ID1(3),2*KL1,2*KL2,2*K,REL)
            RL(KL1+1,KL2+1)=REL
          ELSE
            REL=RL(KL1+1,KL2+1)
          ENDIF
          IF(DABS(REL).GT.EPS) THEN
            S=DBLE((2*KL1+1)*(2*KL2+1)*(2*KS1+1)*(2*KS2+1))
            CA1=HALF*W*RES*REL/DSQRT(S)
          ENDIF
        ENDIF
      ENDIF
C
C     CASES 3231   + + - -        TRANSFORM TO  2133   + - + - 
C
      CALL A1A2W3LSP(IB,IA,IC,IC,KL1,KS1,-HALF,HALF,HALF,-HALF,W)
      IF(DABS(W).LT.EPS) RETURN
      IF(IRS(KS2+1,KS1+1).EQ.0) THEN
        IRS(KS2+1,KS1+1)=1
        CALL RLSP3(3,IB,IA,IC,1,1,2*KS2,2*KS1,2*K,RES)
        RS(KS2+1,KS1+1)=RES
      ELSE
        RES=RS(KS2+1,KS1+1)
      ENDIF
      IF(DABS(RES).GT.EPS)  THEN
        IF(IRL(KL2+1,KL1+1).EQ.0) THEN
          IRL(KL2+1,KL1+1)=1
          CALL 
     :       RLSP3(2,IB,IA,IC,2*ID2(3),2*ID1(3),2*KL2,2*KL1,2*K,REL)
          RL(KL2+1,KL1+1)=REL
        ELSE
          REL=RL(KL2+1,KL1+1)
        ENDIF
        IF(DABS(REL).GT.EPS) THEN
          S=DBLE((2*KL1+1)*(2*KL2+1)*(2*KS1+1)*(2*KS2+1))
          CB1=HALF*W*RES*REL/DSQRT(S)
          IF(MOD(KL1+KL2+KS1+KS2,2).NE.0)CB1=-CB1
        ENDIF
      ENDIF
      RETURN
      END
