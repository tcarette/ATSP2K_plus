*
*     -------------------------------------------------------------
*      N O N R E L A T 5 3 
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE EVALUATED THE CASES - 1423, 4132, 1432, 4123    *
*                                                   ( IREZ = 1),   *
*                                                   ( + + - - ),   *
*     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
*     CONFIGURATIONS:                               N'1 = N1 - 1   *
*                                                   N'2 = N2 + 1   *
*                                                   N'3 = N3 + 1   *
*                                                   N'4 = N4 - 1   *
*     AND    2314, 3241, 2341, 3214                 ( IREZ = 2),   *
*                                                   ( + + - - ),   *
*     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
*     CONFIGURATIONS:                               N'1 = N1 + 1   *
*                                                   N'2 = N2 - 1   *
*                                                   N'3 = N3 - 1   *
*                                                   N'4 = N4 + 1   *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*
      SUBROUTINE NONRELAT53(IA,IB,IC,ID,IREZ)
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
      DIMENSION PMS(2),IATS(2),PML(L2MAX),IATL(L2MAX)
      IF(IHSH.LE.3)RETURN
      IF(.NOT.RECOUPLS0(1,IA,ID,IC,IB,0))RETURN
      IF(.NOT.RECOUPLS0(2,IA,ID,IC,IB,3))RETURN
      IF(.NOT.RECOUPLS0(3,IA,ID,IC,IB,3))RETURN
      IF(IREZ.EQ.1) THEN
        QM1=HALF
        QM2=-HALF
      ELSE
        QM1=-HALF
        QM2=HALF
      ENDIF
      A1=ZERO
      A2=ZERO
      CALL HIBFF(IA,IB,IC,ID,4)
      CALL A1A2A3A4LS(IK1,IK2,IK3,IK4,BK1,BK2,BK3,BK4,ID1,ID2,ID3,ID4,
     :BD1,BD2,BD3,BD4,QM1,QM2,QM2,QM1,ANG)
      IF(DABS(ANG).LT.EPS) RETURN
      PMS(1)=ZERO
      PMS(2)=ZERO
      DO 1 I13=1,2
        JS12=I13-1
        CALL RECOUPLS4(3,IA,IB,IC,ID,1,1,1,1,JS12*2,0,IAT,REC)
        IATS(I13)=IAT
        IF(IATS(I13).NE.0) THEN
          CALL RECOUPLS4(3,IA,IB,IC,ID,1,1,1,1,JS12*2,1,IAT,RECS)
          PMS(I13)=RECS 
        ENDIF
    1 CONTINUE
      IF((IATS(1)+IATS(2)).EQ.0)RETURN
      IP1=ITREXG(LJ(IB),LJ(IA),LJ(IC),LJ(ID),IKK)+1
      IF(IKK.LE.0)RETURN
      IG1=IP1+IKK-1
      LA=IJFUL(IA)
      LB=IJFUL(IB)
      LC=IJFUL(IC)
      LD=IJFUL(ID)
      LIA2=LJ(IA)*2
      LIB2=LJ(IB)*2
      LIC2=LJ(IC)*2
      LID2=LJ(ID)*2
      IAT1=0
      DO 6 I3=IP1,IG1
        PML(I3)=ZERO
        L12=I3-1
        CALL RECOUPLS4(2,IA,IB,IC,ID,LIA2,LIB2,LIC2,LID2,L12*2,0,
     :				IAT,REC)
        IATL(I3)=IAT
        IF(IATL(I3).NE.0) THEN
          CALL RECOUPLS4(2,IA,IB,IC,ID,LIA2,LIB2,LIC2,LID2,L12*2,1,
     :                  IAT,RECL)
          PML(I3)=RECL 
        ENDIF
	IAT1=IAT1+IATL(I3)
    6 CONTINUE
      IF(IAT1.EQ.0)RETURN
      IFAZP=1
      NN=0
      IB1=IB-1
      DO 7 II=IA,IB1
        NN=NOSH1(II)+NN
    7 CONTINUE
      IF((NN/2)*2.EQ.NN)IFAZP=-IFAZP
      NN=0
      ID11=ID-1
      DO 8 II=IC,ID11
        NN=NOSH1(II)+NN
    8 CONTINUE
      IF((NN/2)*2.EQ.NN)IFAZP=-IFAZP
C
C     CASES 1423   + + - -        TRANSFORM TO  1234   + - - +
C           4132                                1234
C                                                    (IREZ = 1)
C     OR
C     CASES 2314   + + - -        TRANSFORM TO  1234   - + + -
C           3241                                1234
C                                                    (IREZ = 2)
C
      IF(IATS(1).NE.0) THEN
        DO 2 I2=IP1,IG1
          KL=I2-1
          IF(CALCULATION(KL)) THEN
            IF((ICOLOM+ISOTOP).EQ.1) THEN
              IF(IATL(I2).NE.0) THEN
                CALL COULOMBLS(LJ(IA),LJ(ID),LJ(IB),LJ(IC),KL,A1)
                IF(DABS(A1).GT.EPS) THEN
                  AB=TWO*A1*ANG*PML(I2)*PMS(1)*DBLE(IFAZP)
     :	                              /DSQRT(DBLE(2*KL+1))
	          IFAZ=IK3(3)+IK4(3)-KL
                  IF(IREZ.EQ.2)IFAZ=IK1(3)+IK2(3)-KL
                  IF((IFAZ/2)*2.NE.IFAZ)AB=-AB
                  IF(DABS(AB).GT.EPS) THEN
                    IF(IREZ.EQ.1)
     :                    CALL SAVENON(3,AB,KL,LA,LD,LB,LC,JA,JB,0)
                    IF(IREZ.EQ.2)
     :                    CALL SAVENON(3,AB,KL,LB,LC,LA,LD,JA,JB,0)
 	          ENDIF
	        ENDIF
	      ENDIF
	    ENDIF
C Orbit 1423
            IF(IORBORB.EQ.1) THEN
              KL2=KL+2
              CALL ORBITORBIT(LJ(IA),LJ(ID),LJ(IB),LJ(IC),KL2,A2)
              IF(DABS(A2).GT.EPS) THEN
                II2=I2+1
                AB=A2*ANG*PML(II2)*PMS(1)*DBLE(IFAZP)
     :                                /DBLE(2*(KL2-1)+1)
                IF(DABS(AB).GT.EPS) THEN
	          IFAZ=IK3(3)+IK4(3)-KL+1
                  IF(IREZ.EQ.2)IFAZ=IK1(3)+IK2(3)-KL+1
                  IF((IFAZ/2)*2.NE.IFAZ)AB=-AB
                  CALL SAVENON(9,AB,KL,LA,LD,LB,LC,JA,JB,0)
                  CALL SAVENON(9,AB,KL,LD,LA,LC,LB,JA,JB,0)
		ENDIF
C                IF(DABS(AB).GT.EPS)
C     :                         WRITE(79,556) AB,KL,LJ(IA),LJ(IB),JA,JB
  556 FORMAT(1X,'1221 O-O','AA=',F17.7,'K=',I3,'LA=',2I3,'JA=',I4,
     :'JB=',I4)
              END IF
            END IF
	  ENDIF
    2   CONTINUE
      ENDIF
C
C     CASES 1432   + + - -        TRANSFORM TO  1234   + - - +
C           4132                                1234
C                                                    (IREZ = 1)
C     OR
C     CASES 2341   + + - -        TRANSFORM TO  1234   - + + -
C           3214                                1234
C                                                    (IREZ = 2)
      IP2=ITREXG(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK)+1
      IF(IKK.LE.0)RETURN
      IG2=IP2+IKK-1
      DO 12 I2=IP2,IG2
        KL=I2-1
        IF(CALCULATION(KL)) THEN
          IF((ICOLOM+ISOTOP).EQ.1) 
     :           CALL COULOMBLS(LJ(IA),LJ(ID),LJ(IC),LJ(IB),KL,A1)
          IF(IORBORB.EQ.1) 
     :           CALL ORBITORBIT(LJ(IA),LJ(ID),LJ(IC),LJ(IB),KL+1,A2)
          IF((DABS(A1)+DABS(A2)).GT.EPS) THEN
            AB=ZERO
            DO 13 I3=IP1,IG1
              L12=I3-1
              IF(IATL(I3).NE.0) THEN
                IF(IXJTIK(LIA2,LIC2,KL*2,LID2,LIB2,L12*2).NE.0) THEN
                  AC=ZERO
                  DO 14 I13=1,2
                    IF(IATS(I13).NE.0) THEN
                      JS12=I13-1
                      AA=PMS(I13)*DSQRT(DBLE(2*JS12+1))
		      AC=AC+AA
	            ENDIF
   14             CONTINUE
                  CALL SIXJ(LIA2,LIC2,KL*2,LID2,LIB2,L12*2,0,SI)
                  AA=AC*PML(I3)*SI*DSQRT(DBLE(2*L12+1))
                  IFAZ=IK3(3)+IK4(3)-L12+1
	          IF(IREZ.EQ.2) IFAZ=IK1(3)+IK2(3)+L12+1
                  IF((IFAZ/2)*2.NE.IFAZ)AA=-AA
                  AB=AB+AA
	        ENDIF
	      ENDIF
   13       CONTINUE
            IF(DABS(AB).GT.EPS) THEN
              AB=ANG*AB*DBLE(IFAZP)
	      ABB=AB
              IF((ICOLOM+ISOTOP).EQ.1) THEN
                IF(DABS(A1).GT.EPS) THEN
                  AB=A1*AB
                  IF(IREZ.EQ.1)
     :                 CALL SAVENON(3,AB,KL,LA,LD,LC,LB,JA,JB,0)
                  IF(IREZ.EQ.2)
     :                 CALL SAVENON(3,AB,KL,LB,LC,LD,LA,JA,JB,0)
 	        ENDIF
 	      ENDIF
C Orbit 1432
              IF(IORBORB.EQ.1) THEN
                IF(DABS(A2).GT.EPS) THEN
	          KLL=KL-1
                  ABB=A2*ABB*HALF/DSQRT(DBLE(2*KL+1))
                  CALL SAVENON(9,ABB,KLL,LA,LD,LC,LB,JA,JB,0)
                  CALL SAVENON(9,ABB,KLL,LD,LA,LB,LC,JA,JB,0)
C                 WRITE(79,655) ABB,KLL,LJ(IA),LJ(IB),JA,JB
  655 FORMAT(1X,'1221 O-O','AA=',F17.7,'K=',I3,'LA=',2I3,'JA=',I4,
     :'JB=',I4)
                END IF
              END IF
 	    ENDIF
          ENDIF
	ENDIF
   12 CONTINUE
      RETURN
      END
