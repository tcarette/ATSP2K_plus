*
*     -------------------------------------------------------------
*      L M A T R I X 2 
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF ONE PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
*                                              N'2 = N2 + 1        *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville           Sebtember 1997   * 
*
      SUBROUTINE LMATRIX2(IA,IB)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      LOGICAL RECOUPLS0
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/DIAGNL/IDIAG,JA,JB
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      IF(IHSH.LE.1)RETURN
      IF(IA.EQ.IB)RETURN
      IF(IA.LT.IB) THEN
         IAA=IA
         IBB=IB
      ELSE
         IAA=IB
         IBB=IA
      ENDIF
      IF(.NOT.RECOUPLS0(1,IAA,IBB,IBB,IBB,0))RETURN
      IF(.NOT.RECOUPLS0(2,IAA,IBB,IBB,IBB,1))RETURN
      IF(.NOT.RECOUPLS0(3,IAA,IBB,IBB,IBB,1))RETURN
      CALL RECOUPLS2(3,IAA,IBB,1,0,IAT,REC)
      IF(IAT.EQ.0)RETURN
      LIA2=LJ(IA)*2
      LIB2=LJ(IB)*2
      CALL RECOUPLS2(2,IAA,IBB,LIB2,0,IAT,REC)
      IF(IAT.EQ.0)RETURN
      CALL RECOUPLS2(3,IAA,IBB,1,1,IAT,RECS)
      CALL RECOUPLS2(2,IAA,IBB,LIB2,1,IAT,RECL)
      CALL HIBFF(IA,IB,IA,IA,2)
      CALL A1A2LS(IK1,IK2,BK1,BK2,ID1,ID2,BD1,BD2,-HALF,HALF,WW)
      IF(DABS(WW).GT.EPS) THEN
        B=WW*RECL*RECS*HALF*DSQRT(DBLE(4*LJ(IA)+2))
        NN=0
        IB1=IBB-1
        DO 1 II=IAA,IB1
          NN=NOSH1(II)+NN
    1   CONTINUE
        IF((NN/2)*2.EQ.NN)B=-B
        IF(DABS(B).GT.EPS) 
     :    CALL SAVENON(4,B,LJ(IA),0,IJFUL(IB),0,IJFUL(IA),JA,JB,0)
      ENDIF
      RETURN
      END
