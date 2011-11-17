*
*     -------------------------------------------------------------
*      W 1
*     -------------------------------------------------------------
*
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
*                                                                  *
*                     N      (K2 K3)  N'      +-                   *
*     ELEMENT:      (l QLS::W      ::l QLS)   -+                   *
*                                             ++                   *
*                                             -- S5(1.47),(1.48),  *
*                                                  (1.49),(1.50).  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE W1(IK,BK,ID,BD,K2,K3,QM1,QM2,W)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DIMENSION IK(7),ID(7),BK(3),BD(3)
      W=ZERO
      IF(ID(3).EQ.3) THEN
	IF(MAX0(IK(4),ID(4)).LT.3) THEN
	  IF(IK(1).LT.300) CALL MES(56)
	  IF(ID(1).LT.300) CALL MES(56)
          IQM2=QM2+QM2+QM2*EPS
          IF((ID(4)+IQM2).GT.2) CALL MES(2)
          CALL W1G(K2,K3,QM1,QM2,IK,BK,ID,BD,W)
          RETURN
        ENDIF
      ELSEIF(ID(3).GT.3) THEN
        IQM2=QM2+QM2+QM2*EPS
        IF((ID(4)+IQM2).GT.2) CALL MES(2)
        CALL W1G(K2,K3,QM1,QM2,IK,BK,ID,BD,W)
        RETURN
      ENDIF
      K=K2+K3
      QQ=QM1+QM2
      IF(ABS(QQ).GE.EPS) THEN
        IF(MOD(K,2).NE.0)RETURN
        IQ=QQ+QQ*TENTH
        IF((IK(4)-ID(4)-2*IQ).NE.0)RETURN
        CALL C1E1SM(BD(1),BD(3),QQ,BK(1),BK(3),A)
        IF(ABS(A).LT.EPS)RETURN
        CALL RWLS(1,K2,K3,IK(3),IK(1),ID(1),AA)
        A=AA*A
        W=A/SQRT(TWO*BK(1)+ONE)
      ELSE
        IF(IK(4).NE.ID(4))RETURN
        IF(K.NE.0) THEN
          K1=1
          IF(MOD(K,2).NE.0)K1=0
          WK1=DBLE(K1)
          CALL CLE0SM(BD(1),BD(3),WK1,BK(1),BK(3),A)
          IF(ABS(A).LT.EPS)RETURN
          CALL RWLS(K1,K2,K3,IK(3),IK(1),ID(1),AA)
          A=AA*A
          W=A/SQRT(FOUR*BK(1)+TWO)
          IF(QM1.GE.EPS)RETURN
          IF(MOD(K,2).NE.0)W=-W
        ELSE
          IF(ID(1).NE.IK(1))RETURN
          IF(QM1.GE.EPS) THEN
            A=-DBLE(ID(4))
          ELSE
            A=DBLE(4*ID(3)+2-ID(4))
          ENDIF
          SS=DBLE((IK(5)+1)*(IK(6)+1))
          S2=DBLE(4*IK(3)+2)
          W=A*SQRT(SS/S2)
        ENDIF
      ENDIF
      RETURN
      END
