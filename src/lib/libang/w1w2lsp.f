*
*     -------------------------------------------------------------
*      W 1 W 2 L S P
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
*     ELEMENTS:                                                    *
*                                                                  *
*        N1       (k1)    N1'                                      * 
*     (nl Q L S::W(11)::nl Q'L'S') *                               *
*        1 1 1 1          1 1 1 1                                  * 
*                                     N2       (k2)    N2'      +- *
*                                * (nl Q L S::W(22)::nl Q'L'S') -+ *
*                                     2 2 2 2          2 2 2 2  ++ *
*                                                               -- *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE W1W2LSP(KL1,KS1,KL2,KS2,QM1,QM2,QM3,QM4,AA)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      COMMON/KAMPAS/ IW1(2,20),IW2(2,20),IWAA(2,2,20,20),
     :RW1(2,20),RW2(2,20),RWAA(2,2,20,20)
      AA=ZERO
      IF(IW1(KS1+1,KL1+1).EQ.0) THEN
        IW1(KS1+1,KL1+1)=1
        IF(IK1(3).GT.3) THEN
          IF(IK1(4).GT.2) CALL MES(35)
          IF(ID1(4).GT.2) CALL MES(35)
        ENDIF
        IQMM1=QM1+QM1+TENTH*QM1
        IQMM2=QM2+QM2+TENTH*QM2
        IF(IK1(4).NE.(ID1(4)+IQMM1+IQMM2)) THEN
          A1=ZERO
        ELSE
          CALL W1(IK1,BK1,ID1,BD1,KL1,KS1,QM1,QM2,A1)
        ENDIF
        RW1(KS1+1,KL1+1)=A1
      ELSE
        A1=RW1(KS1+1,KL1+1)
      ENDIF
      IF(DABS(A1).LT.EPS) RETURN
*
      IF(IW2(KS2+1,KL2+1).EQ.0) THEN
        IW2(KS2+1,KL2+1)=1
        IF(IK2(3).GT.3) THEN
          IF(IK2(4).GT.2) CALL MES(35)
          IF(ID2(4).GT.2) CALL MES(35)
        ENDIF
        IQMM3=QM3+QM3+TENTH*QM3
        IQMM4=QM4+QM4+TENTH*QM4
        IF(IK2(4).NE.(ID2(4)+IQMM3+IQMM4)) THEN
          W=ZERO
        ELSE
          CALL W1(IK2,BK2,ID2,BD2,KL2,KS2,QM3,QM4,W)
        ENDIF
        RW2(KS2+1,KL2+1)=W
      ELSE
        W=RW2(KS2+1,KL2+1)
      ENDIF
      IF(DABS(W).LT.EPS) RETURN
      AA=A1*W
      RETURN
      END
