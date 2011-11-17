*
*     -------------------------------------------------------------
*      W 1 G
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
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE W1G(K1,K2,QM1,QM2,IK,BK,ID,BD,WW)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/RIBOLSF/IMPTLSF(8),IMGTLSF(8),IMPNLSF(8),IMGNLSF(8)
      COMMON/RIBOLS3/IMPTLS3(90),IMGTLS3(90),IMPNLS3(90),IMGNLS3(90)
      DIMENSION IK(7),BK(3),ID(7),BD(3),IBT(7),BT(3)
      WW=ZERO
      IF(IZAS1(ID(7),BD(3),IK(7),BK(3)).EQ.0)RETURN
      ENQP=ZERO
      KK1=K1*2
      KK2=K2*2
      IQ=QM2*TWO+QM2*TENTH
      IF(ID(3).GT.9) RETURN
      IF(ITTK(ID(5),IK(5),KK1).EQ.0)RETURN
      IF(ITTK(ID(6),IK(6),KK2).EQ.0)RETURN
      ITK=IK(1)
      ITD=ID(1)
      IF(ID(3).EQ.3) THEN
        IF(ID(4).GT.2) CALL MES(1)
        IF(IK(4).GT.2) CALL MES(1)
        ITK=ITK-300
        ITD=ITD-300
        ITP1=IMPNLSF(ITK)
        ITP=IMPNLSF(ITD)
        IF(ITP1.NE.ITP)RETURN
        ITG1=IMGNLSF(ITK)
        ITG=IMGNLSF(ITD)
       ELSE
        IF(ID(4).GT.2) CALL MES(1)
        IF(IK(4).GT.2) CALL MES(1)
        ITP1=IMPNLS3(ITK)
        ITP=IMPNLS3(ITD)
        IF(ITP1.NE.ITP)RETURN
        ITG1=IMGNLS3(ITK)
        ITG=IMGNLS3(ITD)
      ENDIF
      IF(ITG1.NE.ITG)RETURN
      IBT(2)=ID(2)
      IBT(3)=ID(3)
      IBT(4)=ID(4)+IQ
      BT(3)=BD(3)+HALF*DBLE(IQ)
      IQM=TWO*DABS(BT(3))+TENTH
      DO 1 IT=ITP,ITG
      CALL RUMT(IT,IBT(3),IBT(7),IBT(6),IBT(5))
      IF(IQM.LE.IBT(7)) THEN
        IF(IXJTIK(IK(3)*2,IK(3)*2,KK1,ID(5),IK(5),IBT(5)).NE.0) THEN
          IF(IXJTIK(1,1,KK2,ID(6),IK(6),IBT(6)).NE.0) THEN
            IBT(1)=IT
            BT(2)=DBLE(IBT(6))/TWO
            BT(1)=DBLE(IBT(7))/TWO
            CALL A1(IK,BK,IBT,BT,QM1,D1)
            IF(DABS(D1).GT.EPS) THEN
              CALL A1(IBT,BT,ID,BD,QM2,W)
              IF(DABS(W).GT.EPS) THEN
                D1=D1*W
                CALL SIXJ(IK(3)*2,IK(3)*2,KK1,ID(5),IK(5),IBT(5),0,SI1)
                CALL SIXJ(1,1,KK2,ID(6),IK(6),IBT(6),0,SI2)
                ENQP=ENQP+D1*SI1*SI2
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    1 CONTINUE
      WW=ENQP*DSQRT(DBLE((KK1+1)*(KK2+1)))
      IE=KK1+IK(6)+ID(6)+KK2+IK(5)+ID(5)
      IF(((IE/4)*4).NE.IE)WW=-WW
      RETURN
      END
