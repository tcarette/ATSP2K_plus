*
*     -------------------------------------------------------------
*      A 1 A 2 W 3 L S P
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
*     ELEMENTS:                                                    *
*                                                                  *
*        N1       (l1,s)   N1'         N2       (l2,s)   N2'       *
*     (nl Q L S:: A(1) ::nl Q'L'S')*(nl Q L S:: A(2) ::nl Q'L'S')* *
*        1 1 1 1           1 1 1 1     2 2 2 2           2 2 2 2   *
*                                                                  *
*         N3        (k1,k2)    N3'                              +- *
*     *(nl Q L S :: W(33) :: nl Q'L'S')                         -+ *
*        3 3 3  3              3 3 3 3                          ++ *
*                                                               -- *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*
      SUBROUTINE A1A2W3LSP(IA,IB,IC,ID,K1,K2,QM1,QM2,QM3,QM4,WW)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      COMMON/TRK2/BD3(3),BD4(3),BK3(3),BK4(3),
     *ID3(7),ID4(7),IK3(7),IK4(7)
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/KAMPAS/ IW1(2,20),IW2(2,20),IWAA(2,2,20,20),
     :RW1(2,20),RW2(2,20),RWAA(2,2,20,20)
      WW=ZERO
      IF(IW1(K2+1,K1+1).EQ.0) THEN
        IW1(K2+1,K1+1)=1
        RW1(K2+1,K1+1)=ZERO
        IF(IK1(3).GT.3) THEN
          IF(IK1(4).GT.2) CALL MES(32)
          IF(ID1(4).GT.2) CALL MES(32)
        ENDIF
        IF(IK2(3).GT.3) THEN
          IF(IK2(4).GT.2) CALL MES(32)
          IF(ID2(4).GT.2) CALL MES(32)
        ENDIF
        IF(IK3(3).GT.3) THEN
          IF(IK3(4).GT.2) CALL MES(32)
          IF(ID3(4).GT.2) CALL MES(32)
        ENDIF
        IQMM1=QM1+QM1+TENTH*QM1
        IF(IK1(4).NE.(ID1(4)+IQMM1))RETURN
        IQMM2=QM2+QM2+TENTH*QM2
        IF(IK2(4).NE.(ID2(4)+IQMM2))RETURN
        IQMM3=QM3+QM3+TENTH*QM3
        IQMM4=QM4+QM4+TENTH*QM4
        IQMM34=IQMM3+IQMM4
        IF(IK3(4).NE.(ID3(4)+IQMM34))RETURN
        CALL C0T5S(BD1(1),BD1(3),QM1,BK1(1),BK1(3),A1)
        IF(DABS(A1).LT.EPS)RETURN
        CALL C0T5S(BD2(1),BD2(3),QM2,BK2(1),BK2(3),C1)
        IF(DABS(C1).LT.EPS)RETURN
        A1=A1*C1
        CALL SLS(IK1(3),IK1(1),IK1(7),IK1(5),IK1(6),ID1(1),ID1(7),
     *                                              ID1(5),ID1(6),S)
        IF(DABS(S).LT.EPS)RETURN
        CALL SLS(IK2(3),IK2(1),IK2(7),IK2(5),IK2(6),ID2(1),ID2(7),
     *                                              ID2(5),ID2(6),C)
        IF(DABS(C).LT.EPS)RETURN
        CALL W1(IK3,BK3,ID3,BD3,K1,K2,QM3,QM4,W)
        IF(DABS(W).LT.EPS)RETURN
        WW=A1*W*S*C/DSQRT(DBLE((IK1(7)+1)*(IK2(7)+1)))
        IFAZP=JFAZE(IA,IB,IC,ID)
	IF(IC.EQ.ID) THEN
          IAAA=MIN0(IA,IB)
          IBBB=MAX0(IA,IB)
	ELSEIF(IA.EQ.ID) THEN
          IAAA=MIN0(IB,IC)
          IBBB=MAX0(IB,IC)
	ELSE
          WRITE(6,'(A)') ' ERROR IN A1A2W3lSP  '
	  STOP
	ENDIF
        NN=1
        IB1=IBBB-1
        DO 1 I=IAAA,IB1
          NN=NOSH1(I)+NN
    1   CONTINUE
        IF(MOD(NN,2).NE.0)IFAZP=-IFAZP
        WW=WW*DBLE(IFAZP)
        RW1(K2+1,K1+1)=WW
      ELSE
        WW=RW1(K2+1,K1+1)
      ENDIF
      RETURN
      END
