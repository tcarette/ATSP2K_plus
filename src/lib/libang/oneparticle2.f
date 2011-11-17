*                                                                  *
*     -------------------------------------------------------------
*      O N E P A R T I C L E 2
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1 -+ 1   *
*                                                  N'2 = N2 +- 1   *
*     Written by G. Gaigalas,                                      *
*     Universite Libre de Bruxelles, Belgium         October 1995  *
*
      SUBROUTINE ONEPARTICLE2(K1,K2,IIA,IIB,XXX)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/DIAGNL/IDIAG,JA,JB
      DIMENSION K(2)
      EXTERNAL XXX
      C1=ZERO
      K(1)=K1
      K(2)=K2
      IA=MIN0(IIA,IIB)
      IB=MAX0(IIA,IIB)
      DO 1 I=1,3
        IF(I.EQ.1) THEN
          CALL RLSP0(I,IA,IB,0,IAT)
        ELSE
          J=I-1
          CALL RLSP0(I,IA,IB,2*K(J),IAT)
        ENDIF
        IF(IAT.EQ.0) RETURN
    1 CONTINUE
      DO 2 I=2,3
        J=I-1
        IF(I.EQ.3) THEN
          CALL RLSP2(I,IA,IB,1,1,2*K(J),0,IAT,REC)
        ELSE
          CALL RLSP2(I,IA,IB,2*LJ(IA),2*LJ(IB),2*K(J),0,IAT,REC)
        ENDIF
        IF(IAT.EQ.0) RETURN
    2 CONTINUE
      CALL HIBFF(IIA,IIB,IIA,IIA,2)
      CALL XXX(IK1(3),IK2(3),INUM,A1)
      IF(DABS(A1).LT.EPS) RETURN
      QM1=HALF
      QM2=-HALF
      IQMM1=QM1+QM1+TENTH*QM1
      IF(IK1(4).NE.(ID1(4)+IQMM1)) RETURN
      IQMM2=QM2+QM2+TENTH*QM2
      IF(IK2(4).NE.(ID2(4)+IQMM2)) RETURN
      CALL C0T5S(BD1(1),BD1(3),QM1,BK1(1),BK1(3),A2)
      IF(DABS(A2).LT.EPS) RETURN
      CALL C0T5S(BD2(1),BD2(3),QM2,BK2(1),BK2(3),A3)
      IF(DABS(A3).LT.EPS) RETURN
      CALL SLS(IK1(3),IK1(1),IK1(7),IK1(5),IK1(6),
     *                       ID1(1),ID1(7),ID1(5),ID1(6),S1)
      IF(DABS(S1).LT.EPS) RETURN
      CALL SLS(IK2(3),IK2(1),IK2(7),IK2(5),IK2(6),
     *                       ID2(1),ID2(7),ID2(5),ID2(6),S2)
      IF(DABS(S2).LT.EPS) RETURN
      LA=IJFUL(IIA)
      LB=IJFUL(IIB)
      A1=S1*S2*A1*A2*A3
      RECLS=1
      DO 3 I=2,3
        J=I-1
        IF(I.EQ.3) THEN
          CALL RLSP2(I,IA,IB,1,1,2*K(J),1,IAT,REC)
        ELSE
          CALL RLSP2(I,IA,IB,2*LJ(IA),2*LJ(IB),2*K(J),1,IAT,REC)
        ENDIF
        RECLS=RECLS*REC
    3 CONTINUE
      C1=A1*RECLS/DSQRT(DBLE((2*K(1)+1)*(2*K(2)+1)
     *                     *(IK1(7)+1)*(IK2(7)+1)))
      JIB=IB-1
      IFAZ=0
      DO 4 I=IA,JIB
        IFAZ=IFAZ+NOSH1(I)
    4 CONTINUE
      IFAZ=IFAZ+1
      IF(MOD(IFAZ,2).NE.0)C1=-C1
      IF(IA.NE.IIA) THEN
        IFAZ=IK1(3)+IK2(3)-K(1)-K(2)
        IF(MOD(IFAZ,2).NE.0)C1=-C1
      ENDIF
      IF(DABS(C1).GT.EPS)CALL SAVENON(INUM,C1,0,0,LA,0,LB,JA,JB,0)
      RETURN
      END
