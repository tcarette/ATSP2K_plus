*
*     -------------------------------------------------------------
*      N O N R E L A T 3 2 
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
*                                              N'2 = N2 + 1        *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*
      SUBROUTINE NONRELAT32(IA,IB,IIA,IIB,IIC,IID)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      LOGICAL RECOUPLS0,CALCULATION
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/DIAGNL/IDIAG,JA,JB
      COMMON /OPERAT/ ICOLOM,ISOTOP,IORBORB
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      IF(IHSH.LE.1)RETURN
      IF(IA.EQ.IB)RETURN
      IF(IA.LT.IB) THEN
         IAA=IA
         IBB=IB
      ELSE
         IAA=IB
         IBB=IA
      ENDIF
      IF(.NOT.RECOUPLS0(1,IAA,IBB,IBB,IBB,0))RETURN
      IF(.NOT.RECOUPLS0(2,IAA,IBB,IBB,IBB,1))RETURN
      IF(.NOT.RECOUPLS0(3,IAA,IBB,IBB,IBB,1))RETURN
      CALL RECOUPLS2(3,IAA,IBB,1,0,IAT,REC)
      IF(IAT.EQ.0)RETURN
      LIA2=LJ(IA)*2
      LIB2=LJ(IB)*2
      CALL RECOUPLS2(2,IAA,IBB,LIA2,0,IAT,REC)
      IF(IAT.EQ.0)RETURN
      CALL RECOUPLS2(3,IAA,IBB,1,1,IAT,RECS)
      CALL RECOUPLS2(2,IAA,IBB,LIA2,1,IAT,RECL)
      QM1=HALF
      QM2=-HALF
      A1=ZERO
      A2=ZERO
      CALL HIBFF(IA,IB,IA,IA,2)
C
C     CASES 2221   + + - -        TRANSFORM TO  1222   - + + -
C           2212                                1222 
C
      LA=IJFUL(IIA)
      LB=IJFUL(IIB)
      LC=IJFUL(IIC)
      LD=IJFUL(IID)
      IP2=ITREXG(LJ(IB),LJ(IB),LJ(IA),LJ(IB),IKK)+1
      IF(IKK.LE.0)RETURN
      IG2=IP2+IKK-1
      DO 2 I2=IP2,IG2
        KL=I2-1
	IF(CALCULATION(KL)) THEN
          IF((ICOLOM+ISOTOP).EQ.1) 
     :           CALL COULOMBLS(LJ(IB),LJ(IB),LJ(IB),LJ(IA),KL,A1)
          IF(IORBORB.EQ.1) 
     :           CALL ORBITORBIT(LJ(IB),LJ(IB),LJ(IB),LJ(IA),KL+1,A2)
          IF((DABS(A1)+DABS(A2)).GT.EPS) THEN
            AB=ZERO
            DO 3 I3=IP2,IG2
              L12=I3-1
              IF(IXJTIK(LIB2,LIA2,KL*2,LIB2,LIB2,L12*2).NE.0) THEN
                AC=ZERO
                DO 4 I13=1,2
                  JS12=I13-1
                  IFAZ=L12+JS12
                  IF((IFAZ/2)*2.EQ.IFAZ) THEN
                    CALL WA1A2LS(IK2,IK1,BK2,BK1,ID2,ID1,BD2,
     *                  BD1,L12,JS12,QM1,QM1,QM2,QM2,AD)
                    AD=AD*DSQRT(DBLE(2*JS12+1))
		    AC=AC+AD
                  ENDIF
    4           CONTINUE
                CALL SIXJ(LIB2,LIA2,KL*2,LIB2,LIB2,L12*2,0,SI)
                AA=AC*SI*DSQRT(DBLE(2*L12+1))
                IFAZ=KL+L12
                IF((IFAZ/2)*2.NE.IFAZ)AA=-AA
                AB=AB+AA
	      ENDIF
    3       CONTINUE
            AB=-AB*RECL*RECS
            IF(DABS(AB).GT.EPS) THEN
	      NN=0
	      IB1=IBB-1
	      DO 5 II=IAA,IB1
	        NN=NOSH1(II)+NN
    5         CONTINUE
	      IF((NN/2)*2.EQ.NN)AB=-AB
	      ABB=AB
              IF((ICOLOM+ISOTOP).EQ.1) THEN
                IF(DABS(A1).GT.EPS) THEN
                  AB=A1*AB
                  CALL SAVENON(3,AB,KL,LA,LB,LC,LD,JA,JB,0)
	        ENDIF
	      ENDIF
C Orbit 2221
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
