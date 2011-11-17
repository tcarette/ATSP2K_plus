*
*     -------------------------------------------------------------
*      W A P 1 G
*     -------------------------------------------------------------
*
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
*                                                                  *
*               N      (k1 K2) N'                    +-            *
*     ELEMENT (l QLS::W     ::l QLS)                 -+            *
*                                                    ++            *
*                                                    -- B17 (2.4)  *
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vilnius,  Lithuania                          December 1993   *
*
      SUBROUTINE WAP1G(K1,K2,K3,BK4,QM1,QM2,QM3,IK,BK,ID,BD,WW)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DIMENSION IK(7),BK(3),ID(7),BD(3),IBT(7),BT(3)
      WW=ZERO
      IF(IZAS1(ID(7),BD(3),IK(7),BK(3)).EQ.0)RETURN
      ENQP=ZERO
      KK1=K1*2
      KK2=K2*2
      KK3=K3*2
      KK4=BK4+BK4+TENTH*BK4
      IQ3=QM3*TWO+QM3*TENTH
      IF(ID(3).GT.9) RETURN
      IF(ITLS3(IK,ID,KK2,KK4,BK,BD,IBT,BT,ITP,ITG,IQ3).EQ.0)RETURN
      IQM=TWO*DABS(BT(3))+TENTH
      DO 1 IT=ITP,ITG
      CALL RUMT(IT,IBT(3),IBT(7),IBT(6),IBT(5))
      IF(IQM.GT.IBT(7))GO TO 1
      IF(IXJTIK(KK1,IK(3)*2,KK2,ID(5),IK(5),IBT(5)).EQ.0)GO TO 1
      IF(IXJTIK(KK3,1,KK4,ID(6),IK(6),IBT(6)).EQ.0)GO TO 1
      IBT(1)=IT
      BT(2)=DBLE(IBT(6))/TWO
      BT(1)=DBLE(IBT(7))/TWO
      CALL A1(IBT,BT,ID,BD,QM3,D1)
      IF(DABS(D1).LT.EPS)GO TO 1
      CALL W1G(K1,K3,QM1,QM2,IK,BK,IBT,BT,W)
      IF(DABS(W).LT.EPS)GO TO 1
      D1=D1*W
      CALL SIXJ(KK1,IK(3)*2,KK2,ID(5),IK(5),IBT(5),0,SI1)
      CALL SIXJ(KK3,1,KK4,ID(6),IK(6),IBT(6),0,SI2)
      ENQP=ENQP+D1*SI1*SI2
    1 CONTINUE
      WW=ENQP*DSQRT(DBLE((KK2+1)*(KK4+1)))
      IE=KK2+IK(6)+ID(6)+KK4+IK(5)+ID(5)
      IF(((IE/4)*4).NE.IE)WW=-WW
      RETURN
      END
