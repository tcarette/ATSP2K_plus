*
*     -------------------------------------------------------------
*      W A 1 A 2 L S P
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
*                                                                  *
*                       N2       (k1,k2) (l2,s) (k3,k4) N2'        *
*     ELEMENT        (nl Q L S::[W(22) * A(2)  ]    ::nl  Q'L'S')* *
*                       2 2 2 2                         2  2 2 2   *
*                                                                  *
*         N1      (l1,s)  N1'                                   +- *
*     *(nl Q L S::A(1)::nl Q'L'S')                              -+ *
*         1 1 1 1         1 1 1 1                               ++ *
*                                                               -- *
*                                                                  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE WA1A2LSP(IAA,IBB,KL1,KS1,KL2,BKKS2,QM1,QM2,QM3,QM4,WW)
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
          IF(IK1(4).GT.2) CALL MES(34)
          IF(ID1(4).GT.2) CALL MES(34)
        ENDIF
        IF(IK2(3).GT.3) THEN
          IF(IK2(4).GT.2) CALL MES(34)
          IF(ID2(4).GT.2) CALL MES(34)
        ENDIF
        IQMM1=QM1+QM1+TENTH*QM1
        IQMM2=QM2+QM2+TENTH*QM2
        IQMM3=QM3+QM3+TENTH*QM3
        IQMM23=IQMM1+IQMM2+IQMM3
        IF(IK2(4).NE.(ID2(4)+IQMM23))RETURN
        IQMM4=QM4+QM4+TENTH*QM4
        IF(IK1(4).NE.(ID1(4)+IQMM4))RETURN
        CALL C0T5S(BD1(1),BD1(3),QM4,BK1(1),BK1(3),A1)
        IF(DABS(A1).LT.EPS)RETURN
        CALL SLS(IK1(3),IK1(1),IK1(7),IK1(5),IK1(6),ID1(1),ID1(7),
     *                                              ID1(5),ID1(6),S)
        IF(DABS(S).LT.EPS)RETURN
        CALL WAP1LS(IK2,BK2,ID2,BD2,KL1,KL2,KS1,BKKS2,QM1,QM2,QM3,WA)
        IF(DABS(WA).LT.EPS)RETURN
        WW=-A1*WA*S/DSQRT(DBLE(IK1(7)+1))
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
