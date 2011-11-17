*
*     -------------------------------------------------------------
*      A 1 A W 2 L S P
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES                           *
*                                        N2      (l2,s)  N2'       *
*     OF FOLLOWING MATRIX ELEMENTS:   (nl Q L S::A(2)::nl Q'L'S')* *
*                                        2 2 2 2         2 2 2 2   *
*                                                                  *
*         N1       (l1,s) (k1,k2) (k3,k4) N1'                   +- *
*     *(nl Q L S::[A(1) * W(11)  ]    ::nl  Q'L'S')             -+ *
*         1 1 1 1                         1  1 1 1              ++ *
*                                                               -- *
*                                                                  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE A1AW2LSP(IAA,IBB,KL1,KS1,KL2,BKKS2,QM1,QM2,QM3,QM4,WW)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      COMMON/KAMPAS/ IW1(2,20),IW2(2,20),IWAA(2,2,20,20),
     :RW1(2,20),RW2(2,20),RWAA(2,2,20,20)
      WW=ZERO
      KKS2=BKKS2+HALF+TENTH
      IF(IWAA(KS1+1,KKS2,KL1+1,KL2+1).EQ.0) THEN
        RWAA(KS1+1,KKS2,KL1+1,KL2+1)=ZERO
        IWAA(KS1+1,KKS2,KL1+1,KL2+1)=1
        IF(IK1(3).GT.3) THEN
          IF(IK1(4).GT.2) CALL MES(33)
          IF(ID1(4).GT.2) CALL MES(33)
        ENDIF
        IF(IK2(3).GT.3) THEN
          IF(IK2(4).GT.2) CALL MES(33)
          IF(ID2(4).GT.2) CALL MES(33)
        ENDIF
        IQMM1=QM1+QM1+TENTH*QM1
        IF(IK2(4).NE.(ID2(4)+IQMM1))RETURN
        IQMM2=QM2+QM2+TENTH*QM2
        IQMM3=QM3+QM3+TENTH*QM3
        IQMM4=QM4+QM4+TENTH*QM4
        IQMM34=IQMM2+IQMM3+IQMM4
        IF(IK1(4).NE.(ID1(4)+IQMM34))RETURN
        CALL C0T5S(BD2(1),BD2(3),QM1,BK2(1),BK2(3),A1)
        IF(DABS(A1).LT.EPS)RETURN
        CALL SLS(IK2(3),IK2(1),IK2(7),IK2(5),IK2(6),ID2(1),ID2(7),
     *                                              ID2(5),ID2(6),S)
        IF(DABS(S).LT.EPS)RETURN
        CALL AWP1LS(IK1,BK1,ID1,BD1,KL1,KL2,KS1,BKKS2,QM2,QM3,QM4,AW)
        IF(DABS(AW).LT.EPS)RETURN
        WW=-A1*AW*S/DSQRT(DBLE(IK2(7)+1))
        NN=1
        IB1=IBB-1
        DO 1 II=IAA,IB1
          NN=NOSH1(II)+NN
    1   CONTINUE
        IF(MOD(NN,2).NE.0)WW=-WW
        RWAA(KS1+1,KKS2,KL1+1,KL2+1)=WW
      ELSE
        WW=RWAA(KS1+1,KKS2,KL1+1,KL2+1)
      ENDIF
      RETURN
      END
