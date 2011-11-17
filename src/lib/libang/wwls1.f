*
*     -------------------------------------------------------------
*      W W L S 1
*     -------------------------------------------------------------
*
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
*                                                                  *
*                N      (k2 K3) (k2 K3) (0 0) N'     +-            *
*     ELEMENT  (l QLS::[W   *  W    ]      ::l QLS)  -+            *
*                                                    ++            *
*                                                    -- B17 (2.4)  *
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vilnius,  Lithuania                             March 1995   *
*
      SUBROUTINE WWLS1(IK,BK,ID,BD,K2,K3,QM1,QM2,QM3,QM4,WW)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DIMENSION IK(7),BK(3),ID(7),BD(3),IBT(7),BT(3)
      WW=ZERO
      IF(ID(5).NE.IK(5))RETURN
      IF(ID(6).NE.IK(6))RETURN
      IF(IZAS1(ID(7),BD(3),IK(7),BK(3)).EQ.0)RETURN
      ENQP=ZERO
      KK2=K2*2
      KK3=K3*2
      IQ3=QM3*TWO+QM3*TENTH
      IQ4=QM4*TWO+QM4*TENTH
      IQ=IQ3+IQ4
      IF(ITLS(IK,ID,0,0,BD,IBT,BT,KK6,KK7,ITP,ITG,IQ).EQ.0)RETURN
      IE1=KK2-IK(5)+KK3-IK(6)
      IQM=TWO*DABS(BT(3))+TENTH
      DO 1 IT=ITP,ITG
      CALL RUMT(IT,IBT(3),IBT(7),IBT(6),IBT(5))
      IF(IQM.LE.IBT(7)) THEN
        IF(IXJTIK(KK2,KK2,0,ID(5),IK(5),IBT(5)).NE.0) THEN
          IF(IXJTIK(KK3,KK3,0,ID(6),IK(6),IBT(6)).NE.0) THEN
            IBT(1)=IT
            BT(2)=DBLE(IBT(6))/TWO
            BT(1)=DBLE(IBT(7))/TWO
            CALL W1(IK,BK,IBT,BT,K2,K3,QM1,QM2,D1)
            IF(DABS(D1).GT.EPS) THEN
              CALL W1(IBT,BT,ID,BD,K2,K3,QM3,QM4,W)
              IF(DABS(W).GT.EPS) THEN
                D1=D1*W
                IF(MOD(IE1+IBT(5)+IBT(6),4).NE.0)D1=-D1
                ENQP=ENQP+D1
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    1 CONTINUE
      WW=ENQP/DSQRT(DBLE((KK2+1)*(IK(5)+1)*(KK3+1)*(IK(6)+1)))
      RETURN
      END
