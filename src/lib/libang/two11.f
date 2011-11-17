*     -------------------------------------------------------------
*      T W O 1 1
*     -------------------------------------------------------------
*                                                                  *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
*                                                                  *
*                                                  1111            *
*                                                                  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE TWO11(KL1,KS1,KL2,KS2,K,IA,C1)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/DIAGNL/IDIAG,JA,JB
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      COMMON/PERMAT/IRS(2,2),IRL(20,20),RS(2,2),RL(20,20)
      C1=ZERO
      W=ZERO
      IF(IXJTIK(2*KL1,2*KL2,2*K,2*ID1(3),2*ID1(3),
     :                                       2*ID1(3)).NE.0) THEN
        IF(IXJTIK(2*KS1,2*KS2,2*K,1,1,1).NE.0) THEN
          CALL W1(IK1,BK1,ID1,BD1,K,K,HALF,-HALF,W)
          IF(DABS(W).GT.EPS) THEN
            CALL SIXJ(2*KL1,2*KL2,2*K,2*ID1(3),
     :                                   2*ID1(3),2*ID1(3),0,SI1)
            CALL SIXJ(1,1,2*K,2*KS1,2*KS2,1,0,SI2)
            W=W*SI1*SI2
          ENDIF
        ENDIF
      ENDIF
      CALL WWPLS1(KL1,KS1,KL2,KS2,K,K,HALF,-HALF,HALF,-HALF,WW)
      WW=WW/DSQRT(DBLE((2*KL1+1)*(2*KL2+1)*(2*KS1+1)*(2*KS2+1)))
      C1=HALF*RS(1,1)*RL(1,1)*(WW-W)
      RETURN
      END
