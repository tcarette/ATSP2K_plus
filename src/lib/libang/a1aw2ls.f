*
*     -------------------------------------------------------------
*      A 1 A W 2 L S
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES                           *
*                                        N1      (l1,s)  N1'       *
*     OF FOLLOWING MATRIX ELEMENTS:   (nl Q L S::A(1)::nl Q'L'S')* *
*                                        1 1 1 1         1 1 1 1   *
*                                                                  *
*         N2       (l2,s) (k1,k2) (l1,s) N2'                    +- *
*     *(nl Q L S::[A(2) * W(22)  ]   ::nl  Q'L'S')              -+ *
*         2 2 2 2                        2  2 2 2               ++ *
*                                                               -- *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*
      SUBROUTINE A1AW2LS(IK1,IK2,BK1,BK2,ID1,ID2,BD1,
     *BD2,K1,K2,QM1,QM2,QM3,QM4,WW)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DIMENSION IK1(7),IK2(7),BK1(3),BK2(3),ID1(7),ID2(7),BD1(3),BD2(3)
      WW=ZERO
      IF(IK1(3).GT.3) THEN
        IF(IK1(4).GT.2) CALL MES(33)
        IF(ID1(4).GT.2) CALL MES(33)
      ENDIF
      IF(IK2(3).GT.3) THEN
        IF(IK2(4).GT.2) CALL MES(33)
        IF(ID2(4).GT.2) CALL MES(33)
      ENDIF
      IQMM1=QM1+QM1+TENTH*QM1
      IF(IK1(4).NE.(ID1(4)+IQMM1))RETURN
      IQMM2=QM2+QM2+TENTH*QM2
      IQMM3=QM3+QM3+TENTH*QM3
      IQMM4=QM4+QM4+TENTH*QM4
      IQMM34=IQMM2+IQMM3+IQMM4
      IF(IK2(4).NE.(ID2(4)+IQMM34))RETURN
      CALL C0T5S(BD1(1),BD1(3),QM1,BK1(1),BK1(3),A1)
      IF(ABS(A1).LT.EPS)RETURN
      CALL SLS(IK1(3),IK1(1),IK1(7),IK1(5),IK1(6),ID1(1),ID1(7),ID1(5)
     *,ID1(6),S)
      IF(ABS(S).LT.EPS)RETURN
      BKK4=HALF
      CALL AWP1LS(IK2,BK2,ID2,BD2,K1,IK1(3),K2,BKK4,QM2,QM3,QM4,AW)
      IF(ABS(AW).LT.EPS)RETURN
      WW=-A1*AW*S/SQRT(DBLE(IK1(7)+1))
      RETURN
      END
