*
*     -------------------------------------------------------------
*      N O N R E L A T 3 3 
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE EVALUATED THE CASES - 2313, 3231, 3213, 2331    *
*                                                   ( + + - - ),   *
*     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
*     CONFIGURATIONS:                               N'1 = N1 - 1   *
*                                                   N'2 = N2 + 1   *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*
      SUBROUTINE NONRELAT33(IA,IB,IC,IREZ,IIA,IIB,IIC,IID)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
*
      PARAMETER(LMAX=10,L2MAX=2*LMAX)
*
      LOGICAL RECOUPLS0,CALCULATION
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/DIAGNL/IDIAG,JA,JB
      COMMON /OPERAT/ ICOLOM,ISOTOP,IORBORB
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      COMMON/TRK2/BD3(3),BD4(3),BK3(3),BK4(3),
     *ID3(7),ID4(7),IK3(7),IK4(7)
      DIMENSION PM(L2MAX,2),IATS(2)
      IF(IHSH.LE.2)RETURN
      CALL EILE(IA,IB,IC,IAA,IBB,ICC)
      IF(.NOT.RECOUPLS0(1,IAA,ICC,IBB,IBB,0))RETURN
      IF(.NOT.RECOUPLS0(2,IAA,ICC,IBB,IBB,2))RETURN
      IF(.NOT.RECOUPLS0(3,IAA,ICC,IBB,IBB,2))RETURN
      LA=IJFUL(IIA)
      LB=IJFUL(IIB)
      LC=IJFUL(IIC)
      LD=IJFUL(IID)
      QM1=HALF
      QM2=-HALF
      A1=ZERO
      A2=ZERO
      CALL HIBFF(IA,IB,IC,IA,3)
      IP1=ITREXG(LJ(IB),LJ(IA),LJ(IC),LJ(IC),IKK)+1
      IF(IKK.LE.0)RETURN
      IG1=IP1+IKK-1
      RECS1=ZERO
      RECS2=ZERO
      CALL RECOUPLS3(3,IB,IA,IC,1,1,0,0,IATS(1),REC)
      IF(IATS(1).NE.0) THEN
        CALL RECOUPLS3(3,IB,IA,IC,1,1,0,1,IAT,RECS1)
      ENDIF
      CALL RECOUPLS3(3,IB,IA,IC,1,1,2,0,IATS(2),REC)
      IF(IATS(2).NE.0) THEN
        CALL RECOUPLS3(3,IB,IA,IC,1,1,2,1,IAT,RECS2)
      ENDIF
      IF((IATS(1)+IATS(2)).EQ.0)RETURN
      LIA2=LJ(IA)*2
      LIB2=LJ(IB)*2
      LIC2=LJ(IC)*2
      DO 4 I4=IP1,IG1
        KL=I4-1
        PM(I4,1)=ZERO
        PM(I4,2)=ZERO
        CALL RECOUPLS3(2,IB,IA,IC,LIB2,LIA2,KL*2,0,IAT,REC)
        IF(IAT.NE.0) THEN
          CALL A1A2W3LS(IK2,IK1,IK3,BK2,BK1,BK3,ID2,ID1,ID3,BD2,
     *                 BD1,BD3,KL,0,QM1,QM2,QM1,QM2,RAG1)
          CALL A1A2W3LS(IK2,IK1,IK3,BK2,BK1,BK3,ID2,ID1,ID3,BD2,
     *                 BD1,BD3,KL,1,QM1,QM2,QM1,QM2,RAG2)
          IF((DABS(RAG1)+DABS(RAG2)).GT.EPS) THEN
            CALL RECOUPLS3(2,IB,IA,IC,LIB2,LIA2,KL*2,1,IAT,REC)
            PM(I4,1)=REC*RECS1*RAG1
            PM(I4,2)=REC*RECS2*RAG2
          ENDIF
        ENDIF
    4 CONTINUE
      IFAZP=JFAZE(IB,IA,IC,IC)
      IF(IA.LT.IB) THEN
        IAA=IA
        IBB=IB
      ELSE
        IAA=IB
        IBB=IA
      ENDIF
      NN=0
      IB1=IBB-1
      DO 16 II=IAA,IB1
        NN=NOSH1(II)+NN
   16 CONTINUE
      IF((NN/2)*2.EQ.NN)IFAZP=-IFAZP
      IF(IREZ.EQ.2)GO TO 5
C
C     CASES 2313   + + - -        TRANSFORM TO  2133   + - + -
C           3231                                2133
C
    6 CONTINUE
      IF(IATS(1).NE.0) THEN
        DO 1 I1=IP1,IG1
          KL=I1-1
	  IF(CALCULATION(KL)) THEN
            IF((ICOLOM+ISOTOP).EQ.1) THEN
              CALL COULOMBLS(LJ(IB),LJ(IC),LJ(IA),LJ(IC),KL,AA)
              AA=AA*PM(I1,1)/DSQRT(DBLE(2*KL+1))
              AA=AA*TWO*DBLE(IFAZP)
              IF(DABS(AA).GT.EPS)
     :                  CALL SAVENON(3,AA,KL,LA,LB,LC,LD,JA,JB,0)
            ENDIF
C Orbit 2313
            IF(IORBORB.EQ.1) THEN
              KL2=KL+2
              CALL ORBITORBIT(LJ(IB),LJ(IC),LJ(IA),LJ(IC),KL2,AA)
              IF(DABS(AA).GT.EPS) THEN
                II1=I1+1
                AA=AA*PM(II1,1)/DBLE(2*(KL2-1)+1)
                AA=AA*DBLE(IFAZP)
                IF(DABS(AA).GT.EPS) THEN
                  CALL SAVENON(9,AA,KL,LA,LB,LC,LD,JA,JB,0)
                  CALL SAVENON(9,AA,KL,LB,LA,LD,LC,JA,JB,0)
                ENDIF
C                IF(DABS(AA).GT.EPS)
C     :                         WRITE(79,556) AA,KL,LJ(IA),LJ(IB),JA,JB
  556 FORMAT(1X,'1221 O-O','AA=',F17.7,'K=',I3,'LA=',2I3,'JA=',I4,
     :'JB=',I4)
              END IF
            END IF
          ENDIF
    1   CONTINUE
      ENDIF
      IF(IREZ.EQ.2)GO TO 7
C
C     CASES 3213   + + - -        TRANSFORM TO  2133   + - + -
C           2331                                2133
C
    5 IP2=ITREXG(LJ(IC),LJ(IA),LJ(IB),LJ(IC),IKK)+1
      IF(IKK.LE.0) GO TO 22
      IG2=IP2+IKK-1
      DO 2 I2=IP2,IG2
        KL=I2-1
	IF(CALCULATION(KL)) THEN
          IF((ICOLOM+ISOTOP).EQ.1) 
     :           CALL COULOMBLS(LJ(IC),LJ(IB),LJ(IA),LJ(IC),KL,A1)
          IF(IORBORB.EQ.1) 
     :           CALL ORBITORBIT(LJ(IC),LJ(IB),LJ(IA),LJ(IC),KL+1,A2)
          IF((DABS(A1)+DABS(A2)).GT.EPS) THEN
            AB=ZERO
            DO 3 I3=IP1,IG1
              L12=I3-1
              IF(IXJTIK(LIA2,LIC2,KL*2,LIC2,LIB2,L12*2).NE.0) THEN
                AC=ZERO
                DO 13 I13=1,2
                  IF(IATS(I13).NE.0) THEN
                    AD=PM(I3,I13)*DSQRT(DBLE(2*(I13-1)+1))
                    IF(I13.EQ.1) AD=-AD
	            AC=AC+AD
                  ENDIF
   13           CONTINUE
                IF(DABS(AC).GT.EPS) THEN
                  CALL SIXJ(LIA2,LIC2,KL*2,LIC2,LIB2,L12*2,0,SI)
                  AA=AC*SI*SQRT(DBLE(2*L12+1))
                  AB=AB+AA
                ENDIF
              ENDIF
    3       CONTINUE
            IF(DABS(AB).GT.EPS) THEN
              AB=AB*DBLE(IFAZP)
	      ABB=AB
              IF((ICOLOM+ISOTOP).EQ.1) THEN
                IF(DABS(A1).GT.EPS) THEN
                  AB=A1*AB
                  CALL SAVENON(3,AB,KL,LA,LB,LD,LC,JA,JB,0)
                ENDIF
              ENDIF
C Orbit 3213
              IF(IORBORB.EQ.1) THEN
                IF(DABS(A2).GT.EPS) THEN
	          KLL=KL-1
                  ABB=A2*ABB*HALF/DSQRT(DBLE(2*KL+1))
                  CALL SAVENON(9,ABB,KLL,LA,LB,LD,LC,JA,JB,0)
                  CALL SAVENON(9,ABB,KLL,LB,LA,LC,LD,JA,JB,0)
C                 WRITE(79,656) ABB,KLL,LJ(IA),LJ(IB),JA,JB
  656 FORMAT(1X,'1221 O-O','AA=',F17.7,'K=',I3,'LA=',2I3,'JA=',I4,
     :'JB=',I4)
                END IF
              END IF
            END IF
          ENDIF
        ENDIF
    2 CONTINUE
   22 IF(IREZ.EQ.2)GO TO 6
    7 CONTINUE
      RETURN
      END
