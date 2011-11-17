*
*     -------------------------------------------------------------
*      W W P L S 1
*     -------------------------------------------------------------
*
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
*                                                                  *
*               N      (k1 K2) (k3 K4) (k5 k6) N'    +-            *
*     ELEMENT (l QLS::[W   *  W    ]        ::l QLS) -+            *
*                                                    ++            *
*                                                    -- B17 (2.4)  *
*                                                                  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE WWPLS1(K1,K2,K3,K4,K5,K6,QM1,QM2,QM3,QM4,WW)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      DIMENSION IBT(7),BT(3)
      WW=ZERO
      IF(IZAS1(ID1(7),BD1(3),IK1(7),BK1(3)).EQ.0)RETURN
      ENQP=ZERO
      KK1=K1*2
      KK2=K2*2
      KK3=K3*2
      KK4=K4*2
      IQ3=QM3*TWO+QM3*TENTH
      IQ4=QM4*TWO+QM4*TENTH
      IQ=IQ3+IQ4
      IF(ITLS(IK1,ID1,K5,K6,BD1,IBT,BT,KK5,KK6,ITP,ITG,IQ).EQ.0)RETURN
      IQM=TWO*DABS(BT(3))+TENTH
      DO 1 IT=ITP,ITG
        CALL RUMT(IT,IBT(3),IBT(7),IBT(6),IBT(5))
        IF(IQM.LE.IBT(7)) THEN
          IF(IXJTIK(KK1,KK3,KK5,ID1(5),IK1(5),IBT(5)).NE.0) THEN
            IF(IXJTIK(KK2,KK4,KK6,ID1(6),IK1(6),IBT(6)).NE.0) THEN
              IBT(1)=IT
              BT(2)=DBLE(IBT(6))/TWO
              BT(1)=DBLE(IBT(7))/TWO
              CALL W1(IK1,BK1,IBT,BT,K1,K2,QM1,QM2,D1)
              IF(DABS(D1).GT.EPS) THEN
                CALL W1(IBT,BT,ID1,BD1,K3,K4,QM3,QM4,W)
                IF(DABS(W).GT.EPS) THEN
                  D1=D1*W
                  CALL SIXJ(KK1,KK3,KK5,ID1(5),IK1(5),IBT(5),0,SI1)
                  CALL SIXJ(KK2,KK4,KK6,ID1(6),IK1(6),IBT(6),0,SI2)
                  ENQP=ENQP+D1*SI1*SI2
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
    1 CONTINUE
      WW=ENQP*DSQRT(DBLE((KK5+1)*(KK6+1)))
      IF(MOD(KK5+IK1(6)+ID1(6)+KK6+IK1(5)+ID1(5),4).NE.0) WW=-WW
      RETURN
      END
