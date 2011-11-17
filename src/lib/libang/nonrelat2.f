*
*     -------------------------------------------------------------
*      N O N R E L A T 2 
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 2        *
*                                              N'2 = N2 + 2        *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*
      SUBROUTINE NONRELAT2(IA,IB)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      LOGICAL IATT,IATTT,RECOUPLS0,CALCULATION
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/DIAGNL/IDIAG,JA,JB
      COMMON /OPERAT/ ICOLOM,ISOTOP,IORBORB
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      DIMENSION PMS(2),IATT(2) 
      IF(IHSH.LE.1)RETURN
      IF(IA.EQ.IB)RETURN
        IF(ISOTOP.EQ.1) THEN
	  IF(ITTK(LJ(IA),LJ(IB),1).EQ.0)RETURN
        ENDIF
      IF(IA.LT.IB) THEN
         IAA=IA
         IBB=IB
      ELSE
         IAA=IB
         IBB=IA
      ENDIF
      IF(.NOT.RECOUPLS0(1,IAA,IBB,IBB,IBB,0))RETURN
      IF(.NOT.RECOUPLS0(2,IAA,IBB,IBB,IBB,1))RETURN
      IATT(1)=RECOUPLS0(3,IAA,IBB,IBB,IBB,0)
      IATT(2)=RECOUPLS0(3,IAA,IBB,IBB,IBB,1)
      IATTT=.TRUE.
      IF(IATT(1).OR.IATT(2))IATTT=.FALSE.
      IF(IATTT)RETURN
      QM1=HALF
      QM2=-HALF
      A1=ZERO
      A2=ZERO
      CALL HIBFF(IA,IB,IA,IA,2)
C
C     THE CASE 1122   + + - -
C
      PMS(1)=ZERO
      PMS(2)=ZERO
      DO 1 I13=1,2
	IF(IATT(I13)) THEN
          JS12=I13-1
          CALL RECOUPLS2(3,IAA,IBB,JS12*2,0,IAT,REC)
          IF(IAT.NE.0) THEN
            CALL RECOUPLS2(3,IAA,IBB,JS12*2,1,IAT,RECS)
            PMS(I13)=RECS 
          ENDIF
	ENDIF
    1 CONTINUE
      LA=IJFUL(IA)
      LB=IJFUL(IB)
      LIA2=LJ(IA)*2
      LIB2=LJ(IB)*2
      IG1=MIN(LIA2,LIB2)+1
      IP2=IABS(LJ(IB)-LJ(IA))+1
      IG2=LJ(IB)+LJ(IA)+1
      IGAL=2
      IF(IORBORB.EQ.1) IGAL=1
      DO 2 I2=IP2,IG2,IGAL
        KL=I2-1
	IF(CALCULATION(KL)) THEN
          IF((ICOLOM+ISOTOP).EQ.1) 
     :           CALL COULOMBLS(LJ(IA),LJ(IA),LJ(IB),LJ(IB),KL,A1)
          IF(IORBORB.EQ.1) 
     :           CALL ORBITORBIT(LJ(IA),LJ(IA),LJ(IB),LJ(IB),KL+1,A2)
          IF((DABS(A1)+DABS(A2)).GT.EPS) THEN
            A1=-HALF*A1
            AB=ZERO
            DO 3 I3=1,IG1
              L12=I3-1
              CALL RECOUPLS2(2,IAA,IBB,L12*2,0,IAT,REC)
              IF(IAT.NE.0) THEN
                IF(IXJTIK(LIA2,LIB2,KL*2,LIB2,LIA2,L12*2).NE.0) THEN
	          AC=ZERO
                  DO 4 I13=1,2
	            IF(IATT(I13)) THEN
                      JS12=I13-1
                      CALL W1W2LS(L12,JS12,L12,JS12,QM1,QM1,QM2,QM2,AA)
                      AA=AA*PMS(I13)*DSQRT(DBLE(2*JS12+1))
		      AC=AC+AA
	            ENDIF
    4             CONTINUE
                  CALL RECOUPLS2(2,IAA,IBB,L12*2,1,IAT,RECL)
                  CALL SIXJ(LIA2,LIB2,KL*2,LIB2,LIA2,L12*2,0,SI)
                  AA=AC*RECL*SI*DSQRT(DBLE(2*L12+1))
                  IFAZ=IK1(3)+IK2(3)+KL+L12
                  IF((IFAZ/2)*2.NE.IFAZ)AA=-AA
                  AB=AB+AA
	        ENDIF
	      ENDIF
    3       CONTINUE
            IF(DABS(AB).GT.EPS) THEN
	      ABB=AB
              IF((ICOLOM+ISOTOP).EQ.1) THEN
                AB=A1*AB
                IF(DABS(AB).GT.EPS)
     :	            CALL SAVENON(3,AB,KL,LA,LA,LB,LB,JA,JB,0)
              ENDIF
C Orbit 1122
              IF(IORBORB.EQ.1) THEN
                IF(DABS(A2).GT.EPS) THEN
	          KLL=KL-1
                  ABB=-A2*ABB*HALF/DSQRT(DBLE(2*KL+1))
                  CALL SAVENON(9,ABB,KLL,LA,LA,LB,LB,JA,JB,0)
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
