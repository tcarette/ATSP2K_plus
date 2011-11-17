*
*     -------------------------------------------------------------
*      A W P 1 L S
*     -------------------------------------------------------------
*
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
*                                                                  *
*                N      (ls)  (k1 k3) (k2 k4)  N'     +-           *
*     ELEMENT: (l QLS ::[A  * W      ]      ::l  QLS) -+           *
*                                                     ++           *
*                                                     -- B17 (2.3) *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE AWP1LS(IK,BK,ID,BD,K1,K2,K3,BK4,QM1,QM2,QM3,AW)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DIMENSION IK(7),BK(3),ID(7),BD(3),IBT(7),BT(3)
      AW=ZERO
      IF(ID(3).EQ.3) THEN
	IF(MAX0(IK(4),ID(4)).LT.3) THEN
	  IF(IK(1).LT.300) CALL MES(54)
	  IF(ID(1).LT.300) CALL MES(54)
          CALL AWP1G(K1,K2,K3,BK4,QM1,QM2,QM3,IK,BK,ID,BD,AW)
          RETURN
        ENDIF
      ELSEIF(ID(3).GT.3) THEN
        CALL AWP1G(K1,K2,K3,BK4,QM1,QM2,QM3,IK,BK,ID,BD,AW)
        RETURN
      ENDIF
      IF(IZAS1(ID(7),BD(3),IK(7),BK(3)).EQ.0)RETURN
      ENQP=ZERO
      IQ2=QM2*TWO+QM2*TENTH
      IQ3=QM3*TWO+QM3*TENTH
      IQ=IQ2+IQ3
      KK1=K1*2
      KK2=K2*2
      KK3=K3*2
      KK4=BK4+BK4+TENTH*BK4
      IF(ITLS2(IK,ID,KK2,KK4,BD,IBT,BT,ITP,ITG,IQ).EQ.0)RETURN
      IQM=TWO*DABS(BT(3))+TENTH
      DO 1 IT=ITP,ITG
        CALL RUMT(IT,IBT(3),IBT(7),IBT(6),IBT(5))
        IF(IQM.GT.IBT(7))GO TO 1
        IF(IXJTIK(2*IK(3),KK1,KK2,ID(5),IK(5),IBT(5)).EQ.0)GO TO 1
        IF(IXJTIK(1,KK3,KK4,ID(6),IK(6),IBT(6)).EQ.0)GO TO 1
        IBT(1)=IT
        BT(2)=DBLE(IBT(6))/TWO
        BT(1)=DBLE(IBT(7))/TWO
        CALL C0T5S(BT(1),BT(3),QM1,BK(1),BK(3),D1)
        IF(DABS(D1).LT.EPS)GO TO 1
        CALL SLS(IK(3),IK(1),IK(7),IK(5),IK(6),IBT(1),IBT(7),IBT(5),
     *  IBT(6),S)
        IF(DABS(S).LT.EPS)GO TO 1
        CALL W1(IBT,BT,ID,BD,K1,K3,QM2,QM3,W)
        IF(DABS(W).LT.EPS)GO TO 1
        D1=D1*W*S
        CALL SIXJ(2*IK(3),KK1,KK2,ID(5),IK(5),IBT(5),0,SI1)
        CALL SIXJ(1,KK3,KK4,ID(6),IK(6),IBT(6),0,SI2)
        D1=D1*SI1*SI2/DSQRT(DBLE(IK(7)+1))
        ENQP=ENQP+D1
    1 CONTINUE
      AW=ENQP*DSQRT(DBLE((KK2+1)*(KK4+1)))
      IE=KK2+IK(6)+ID(6)+2+KK4+IK(5)+ID(5)
      IF(((IE/4)*4).NE.IE)AW=-AW
      RETURN
      END
