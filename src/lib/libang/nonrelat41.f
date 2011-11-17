*
*     -------------------------------------------------------------
*      N O N R E L A T 4 1 
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 + 1        *
*                                              N'2 = N2 + 1        *
*                                              N'3 = N3 - 2,       *
*     WHEN IREZ = 1   ........................
*                                              N'1 = N1 - 1        *
*                                              N'2 = N2 - 1        *
*                                              N'3 = N3 + 2,       *
*     WHEN IREZ = 2   ........................
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*
      SUBROUTINE NONRELAT41(IA,IB,IC,IREZ,IIA,IIB,IIC,IID)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
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
      DIMENSION PMS(2),IATS(2)
      IF(IHSH.LE.2)RETURN
      CALL EILE(IA,IB,IC,IAA,IBB,ICC)
      IF(.NOT.RECOUPLS0(1,IAA,ICC,IBB,IBB,0))RETURN
      IF(.NOT.RECOUPLS0(2,IAA,ICC,IBB,IBB,2))RETURN
      IF(.NOT.RECOUPLS0(3,IAA,ICC,IBB,IBB,2))RETURN
      IF(IREZ.EQ.1) THEN
        QM1=-HALF
        QM2=HALF
      ELSE
        QM1=HALF
        QM2=-HALF
      ENDIF
      A1=ZERO
      A2=ZERO
      CALL HIBFF(IA,IB,IC,IA,3)
C
C     CASES 3312   + + - -        TRANSFORM TO  1233   - - + +
C           3321                                1233
C                                                    (IREZ = 1)
C     OR
C     CASES 1233   + + - -        TRANSFORM TO  1233   + + - -
C           2133                                1233
C                                                    (IREZ = 2)
      PMS(1)=ZERO
      PMS(2)=ZERO
      DO 1 I13=1,2
        JS12=I13-1
        CALL RECOUPLS3(3,IA,IB,IC,1,1,JS12*2,0,IAT,REC)
        IATS(I13)=IAT
        IF(IATS(I13).NE.0) THEN
          CALL RECOUPLS3(3,IA,IB,IC,1,1,JS12*2,1,IAT,RECS)
          PMS(I13)=RECS 
        ENDIF
    1 CONTINUE
      IF((IATS(1)+IATS(2)).EQ.0)RETURN
      LA=IJFUL(IIA)
      LB=IJFUL(IIB)
      LC=IJFUL(IIC)
      LD=IJFUL(IID)
      LIA2=LJ(IA)*2
      LIB2=LJ(IB)*2
      LIC2=LJ(IC)*2
      IP1=ITREXG(LJ(IB),LJ(IA),LJ(IC),LJ(IC),IKK)+1
      IF(IKK.LE.0)RETURN
      IG1=IP1+IKK-1
      IP2=ITREXG(LJ(IC),LJ(IA),LJ(IB),LJ(IC),IKK)+1
      IF(IKK.LE.0)RETURN
      IG2=IP2+IKK-1
      IFAZP=JFAZE(IC,IA,IB,IC)
      IF(IA.LT.IB) THEN
         IAA=IA
         IBB=IB
      ELSE
         IAA=IB
         IBB=IA
      ENDIF
      NN=0
      IB1=IBB-1
      DO 5 II=IAA,IB1
        NN=NOSH1(II)+NN
    5 CONTINUE
      IF((NN/2)*2.EQ.NN)IFAZP=-IFAZP
      DO 2 I2=IP2,IG2
        KL=I2-1
	IF(CALCULATION(KL)) THEN
          IF((ICOLOM+ISOTOP).EQ.1) 
     :           CALL COULOMBLS(LJ(IA),LJ(IB),LJ(IC),LJ(IC),KL,A1)
          IF(IORBORB.EQ.1) 
     :           CALL ORBITORBIT(LJ(IA),LJ(IB),LJ(IC),LJ(IC),KL+1,A2)
          IF((DABS(A1)+DABS(A2)).GT.EPS) THEN
            AB=ZERO
            DO 3 I3=IP1,IG1
              L12=I3-1
              CALL RECOUPLS3(2,IA,IB,IC,LIA2,LIB2,L12*2,0,IAT,REC)
              IF(IAT.NE.0) THEN
                IF(IXJTIK(LIC2,LIB2,KL*2,LIA2,LIC2,L12*2).NE.0) THEN
                  AC=ZERO
                  DO 4 I13=1,2
		    IF(IATS(I13).NE.0) THEN
                      JS12=I13-1
		      IFAZ=L12+JS12
	              IF((IFAZ/2)*2.EQ.IFAZ) THEN
                        CALL A1A2W3LS(IK1,IK2,IK3,BK1,BK2,BK3,
     :          ID1,ID2,ID3,BD1,BD2,BD3,L12,JS12,QM1,QM1,QM2,QM2,AA)
                        AA=AA*PMS(I13)*DSQRT(DBLE(2*JS12+1))
		        AC=AC+AA
	              ENDIF
	            ENDIF
    4             CONTINUE
                  CALL RECOUPLS3(2,IA,IB,IC,LIA2,LIB2,L12*2,1,IAT,RECL)
                  CALL SIXJ(LIC2,LIB2,KL*2,LIA2,LIC2,L12*2,0,SI)
                  AA=AC*RECL*SI*DSQRT(DBLE(2*L12+1))
                  IFAZ=IK1(3)+IK3(3)+KL+L12
	          IF(IREZ.EQ.2) IFAZ=IK2(3)+IK3(3)+KL+L12
                  IF((IFAZ/2)*2.NE.IFAZ)AA=-AA
                  AB=AB+AA
	        ENDIF
	      ENDIF
    3       CONTINUE
            IF(DABS(AB).GT.EPS) THEN
              AB=-AB*DBLE(IFAZP)
	      ABB=AB
              IF((ICOLOM+ISOTOP).EQ.1) THEN
                IF(DABS(A1).GT.EPS) THEN
                  AB=A1*AB
                  CALL SAVENON(3,AB,KL,LA,LB,LC,LD,JA,JB,0)
                ENDIF
              ENDIF
C Orbit 3312
              IF(IORBORB.EQ.1) THEN
                IF(DABS(A2).GT.EPS) THEN
	          KLL=KL-1
                  ABB=A2*ABB*HALF/DSQRT(DBLE(2*KL+1))
                  CALL SAVENON(9,ABB,KLL,LA,LB,LC,LD,JA,JB,0)
                  CALL SAVENON(9,ABB,KLL,LB,LA,LD,LC,JA,JB,0)
C                 WRITE(79,656) ABB,KLL,LJ(IA),LJ(IB),JA,JB
  656 FORMAT(1X,'1221 O-O','AA=',F17.7,'K=',I3,'LA=',2I3,'JA=',I4,
     :'JB=',I4)
                END IF
              END IF
            ENDIF
          ENDIF
	ENDIF
    2 CONTINUE
      RETURN
      END
