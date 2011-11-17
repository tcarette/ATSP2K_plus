*                                                                  *
*     -------------------------------------------------------------
*      T W O 1 2
*     -------------------------------------------------------------
*                                                                  *
*                                                                  *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
*                                                  N'2 = N2        *
*                                                                  *
*           1212   + + - -        TRANSFORM TO  1122   + - + -     *
*           2121                                1122               *
*                                                                  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE TWO12(KL1,KS1,KL2,KS2,K,IA,IB,C1)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON /CASEOP/ IOCASE
      COMMON/PERMAT/IRS(2,2),IRL(20,20),RS(2,2),RL(20,20)
      C1=ZERO
      IF(IRS(KS1+1,KS2+1).EQ.0) THEN
	IRS(KS1+1,KS2+1)=1
        CALL RLSP2(3,IA,IB,2*KS1,2*KS2,2*K,0,IAT,RECS)
        IF(IAT.EQ.0) THEN
          RECS=ZERO
        ELSE
          CALL RLSP2(3,IA,IB,2*KS1,2*KS2,2*K,1,IAT,RECS)
        ENDIF
	RS(KS1+1,KS2+1)=RECS
      ELSE
	RECS=RS(KS1+1,KS2+1)
      ENDIF
      IF(DABS(RECS).LT.EPS) RETURN
      IF(IRL(KL1+1,KL2+1).EQ.0) THEN
	IRL(KL1+1,KL2+1)=1
        CALL RLSP2(2,IA,IB,2*KL1,2*KL2,2*K,0,IAT,RECL)
        IF(IAT.EQ.0) THEN
          RECL=ZERO
        ELSE
          CALL RLSP2(2,IA,IB,2*KL1,2*KL2,2*K,1,IAT,RECL)
        ENDIF
	RL(KL1+1,KL2+1)=RECL
      ELSE
	RECL=RL(KL1+1,KL2+1)
      ENDIF
      IF(DABS(RECL).LT.EPS) RETURN
      CALL W1W2LSP(KL1,KS1,KL2,KS2,HALF,-HALF,HALF,-HALF,W)
      IF(DABS(W).LT.EPS) RETURN
      W=W/DSQRT(DBLE((2*KL1+1)*(2*KL2+1)*(2*KS1+1)*(2*KS2+1)))
      RECLS=RECL*RECS
      C1=RECLS*W*HALF
      RETURN
      END
