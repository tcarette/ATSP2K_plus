*
*     -------------------------------------------------------------
*      N O N R E L A T 1 
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :           N'1 = N1        *
*                                                  N'2 = N2        *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*                                                                  * 
*     Universite Libre de Bruxelles, Brussels                      * 
*                                                  December 1995   * 
*
      SUBROUTINE NONRELAT1(IA,IB,IIRE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
*
      PARAMETER(LMAX=10,L2MAX=2*LMAX)
*
      LOGICAL IATT,IATTT,RECOUPLS0,CALCULATION
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/DIAGNL/IDIAG,JA,JB
      COMMON /OPERAT/ ICOLOM,ISOTOP,IORBORB
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      DIMENSION PMG(L2MAX,2),RAG(L2MAX,2),IATT(2)
      LIA2=LJ(IA)*2
      LA=IJFUL(IA)
      QM1=HALF
      QM2=-HALF
      A1=ZERO
      A2=ZERO
      IF(IA.EQ.IB) THEN
C
C     THE CASE 1111   + + - -
C
        IF(JA.NE.JB)THEN
          DO 1 IN=1,3
            IF(.NOT.RECOUPLS0(IN,IA,IA,IA,IA,0))RETURN
    1     CONTINUE
        END IF
        CALL HIBFF(IA,IA,IA,IA,1)
        CALL W1(IK1,BK1,ID1,BD1,0,0,QM1,QM2,A)
        RECOPL=ONE/DSQRT(DBLE(IK1(6)+1)*(IK1(5)+1))
        IF(ABS(A).GT.EPS) THEN
          B=A*RECOPL*HALF*DSQRT(DBLE(4*LJ(IA)+2))
CBAY            IF(DABS(B).GT.EPS) CALL SAVENON(4,B,0,0,LA,0,LA,JA,JB,0)
          IF (DABS(B).GT.EPS) THEN
            if(ja.ne.jb) then
              write(*,'(A,I5,A,I5)')
     :            'Configurations',ja,'and',jb,'are identical'
              STOP
            ELSE
              CALL SAVENON(4,B,0,0,LA,0,LA,JA,JB,0)
            END IF
          END IF
          A=A/DSQRT(DBLE(4*LJ(IA)+2))
        END IF
        IF(IIRE.EQ.0)RETURN
        IF((ICOLOM+IORBORB).NE.0) THEN
          IG2=LIA2+1
          DO 2 I2=1,IG2,2
            KL=I2-1
            IF(ICOLOM.EQ.1) THEN
CGGf                                with    f-shell   ***  beginning
              IF(JA.NE.JB)THEN
CGGf                                with    f-shell   ***  e n d
                CALL COULOMBLS(LJ(IA),LJ(IA),LJ(IA),LJ(IA),KL,A1)
                IF(DABS(A1).GT.EPS) THEN
                  CALL WWLS1(IK1,BK1,ID1,BD1,KL,0,QM1,QM2,QM1,QM2,AA)
                  AA=(AA/DSQRT(DBLE(2*KL+1)))+A
	          AA=AA*A1*RECOPL
                  IF(DABS(AA).GT.EPS)
     :                      CALL SAVENON(1,AA,KL,0,LA,0,LA,JA,JB,0)
                END IF
              ELSEIF(LJ(IA).GT.3) THEN
                CALL COULOMBLS(LJ(IA),LJ(IA),LJ(IA),LJ(IA),KL,A1)
                IF(DABS(A1).GT.EPS) THEN
                  CALL WWLS1(IK1,BK1,ID1,BD1,KL,0,QM1,QM2,QM1,QM2,AA)
                  AA=(AA/DSQRT(DBLE(2*KL+1)))+A
	          AA=AA*A1*RECOPL
                  IF(DABS(AA).GT.EPS)
     :                      CALL SAVENON(1,AA,KL,0,LA,0,LA,JA,JB,0)
                END IF
CGGf                                with    f-shell   ***  beginning
              ELSE
                CALL AVERA(ID1(3),ID1(1),KL,ID1(4),AA)
                IF(DABS(AA).GT.EPS)
     :                    CALL SAVENON(1,AA,KL,0,LA,0,LA,JA,JB,0)
CGGf                                with    f-shell   ***  e n d
              END IF
            END IF
C Orbit 1111
            IF(IORBORB.EQ.1) THEN
              KL2=KL+2
              KL1=KL+1
              CALL ORBITORBIT(LJ(IA),LJ(IA),LJ(IA),LJ(IA),KL2,A1)
              IF(DABS(A1).GT.EPS) THEN
                CALL WWLS1(IK1,BK1,ID1,BD1,KL1,0,QM1,QM2,QM1,QM2,AA)
                AA=(AA/DBLE(2*KL1+1))-A/DSQRT(DBLE(2*KL1+1))
	        AA=AA*A1*RECOPL
                IF(DABS(AA).GT.EPS)
     :             CALL SAVENON(9,AA,KL,LA,LA,LA,LA,JA,JB,0)
C                IF(DABS(AA).GT.EPS)
C     :                               WRITE(77,555) AA,KL,LA,JA,JB
  555 FORMAT(1X,'1111 O-O','AA=',F17.7,'K=',I3,'LA=',I3,'J=',2I4)
              END IF
            END IF
    2     CONTINUE
          RETURN
        END IF
      END IF
      IF(IIRE.EQ.0)RETURN
      IF(IHSH.LE.1)RETURN
      IF(ISOTOP.EQ.1) THEN
        IF(ITTK(LJ(IA),LJ(IB),1).EQ.0)RETURN
      ENDIF
      IATT(1)=.TRUE.
      IATT(2)=.TRUE.
      IF(JA.NE.JB)THEN
	 IF(.NOT.RECOUPLS0(1,IA,IB,IB,IB,0))RETURN
	 IF(.NOT.RECOUPLS0(2,IA,IB,IB,IB,1))RETURN
         IATT(1)=RECOUPLS0(3,IA,IB,IB,IB,0)
         IATT(2)=RECOUPLS0(3,IA,IB,IB,IB,1)
	 IATTT=.TRUE.
         IF(IATT(1).OR.IATT(2))IATTT=.FALSE.
         IF(IATTT)RETURN
      END IF
      CALL HIBFF(IA,IB,IA,IA,2)
      LIB2=LJ(IB)*2
      LB=IJFUL(IB)
      IG1=MIN(LIA2,LIB2)+1
      DO 3 I4=1,IG1
        KL=I4-1
        RAG(I4,1)=ZERO
        PMG(I4,1)=ZERO
        RAG(I4,2)=ZERO
        PMG(I4,2)=ZERO
        CALL RECOUPLS2(2,IA,IB,KL*2,0,IAT,RE)
        IF(IAT.NE.0) THEN
          CALL RECOUPLS2(2,IA,IB,KL*2,1,IAT,RE)
          IF(IATT(1)) THEN
            CALL RECOUPLS2(3,IA,IB,0,0,IAT,REP)
            IF(IAT.NE.0) THEN
              CALL W1W2LS(KL,0,KL,0,QM1,QM2,QM1,QM2,RA)
              RAG(I4,1)=RA
              PMG(I4,1)=RE/DSQRT(DBLE((IK1(6)+1)*(IK2(6)+1)))
	    END IF
	  END IF
          IF(IATT(2)) THEN
            CALL RECOUPLS2(3,IA,IB,2,0,IAT,REP)
            IF(IAT.NE.0) THEN
              CALL W1W2LS(KL,1,KL,1,QM1,QM2,QM1,QM2,RA)
              IF(DABS(RA).GT.EPS) THEN
                RAG(I4,2)=RA
                CALL RECOUPLS2(3,IA,IB,2,1,IAT,REP)
                PMG(I4,2)=RE*REP
	      END IF
	    END IF
	  END IF
	END IF
    3 CONTINUE
C 
C     CASES 1212   + + - -        TRANSFORM TO  1122   + - + -
C           2121                                1122
C
      IF(IATT(1)) THEN
        IF((ICOLOM+IORBORB).NE.0) THEN
          DO 4 I1=1,IG1,2
            KL=I1-1
            IF(ICOLOM.EQ.1) THEN
              CALL COULOMBLS(LJ(IA),LJ(IB),LJ(IA),LJ(IB),KL,AA)
              IF(DABS(AA).GE.EPS) THEN
                AA=AA*PMG(I1,1)
                AA=AA*RAG(I1,1)
                AA=AA*TWO/DSQRT(DBLE(2*KL+1))
                IF(DABS(AA).GT.EPS)
     :              CALL SAVENON(1,AA,KL,0,LA,0,LB,JA,JB,0)
              END IF
            END IF
C Orbit 1212
            IF(IORBORB.EQ.1) THEN
              KL2=KL+2
              CALL ORBITORBIT(LJ(IA),LJ(IB),LJ(IA),LJ(IB),KL2,AA)
              IF(DABS(AA).GT.EPS) THEN
                II1=I1+1
                AA=AA*PMG(II1,1)
                AA=AA*RAG(II1,1)
                AA=AA/DBLE(2*(KL2-1)+1)
                IF(DABS(AA).GT.EPS) THEN
                  CALL SAVENON(9,AA,KL,LA,LB,LA,LB,JA,JB,0)
                  CALL SAVENON(9,AA,KL,LB,LA,LB,LA,JA,JB,0)
		ENDIF
C                IF(DABS(AA).GT.EPS)
C     :                         WRITE(78,556) AA,KL,LJ(IA),LJ(IB),JA,JB
  556 FORMAT(1X,'1212 O-O','AA=',F17.7,'K=',I3,'LA=',2I3,'JA=',I4,
     :'JB=',I4)
              END IF
            END IF
    4     CONTINUE
        END IF
      END IF
C
C     CASES 1221   + + - -        TRANSFORM TO  1122   + - + -
C           2112                                1122
C
      IP2=IABS(LJ(IB)-LJ(IA))+1
      IG2=LJ(IB)+LJ(IA)+1
      IGAL=2
      IF(IORBORB.EQ.1) IGAL=1
      DO 5 I2=IP2,IG2,IGAL
        KL=I2-1
	IF(CALCULATION(KL)) THEN
          IF((ICOLOM+ISOTOP).EQ.1) 
     :           CALL COULOMBLS(LJ(IA),LJ(IB),LJ(IB),LJ(IA),KL,A1)
          IF(IORBORB.EQ.1) 
     :           CALL ORBITORBIT(LJ(IA),LJ(IB),LJ(IB),LJ(IA),KL+1,A2)
          IF((DABS(A1)+DABS(A2)).GT.EPS) THEN
            AB=ZERO
            DO 6 I3=1,IG1
              L12=I3-1
              AC=ZERO
              IF(IXJTIK(LIA2,LIB2,KL*2,LIB2,LIA2,L12*2).NE.0)THEN
                DO 7 I13=1,2
                  IF(IATT(I13)) THEN
                    AA=PMG(I3,I13)
                    AA=AA*RAG(I3,I13)
                    AA=AA*DSQRT(DBLE(2*(I13-1)+1))
                    IF(I13.EQ.1)AA=-AA
                    AC=AC+AA
                  END IF
    7           CONTINUE
                IF(DABS(AC).GT.EPS) THEN 
                  CALL SIXJ(LIA2,LIB2,KL*2,LIB2,LIA2,L12*2,0,SI)
                  AA=AC*SI*DSQRT(DBLE(2*L12+1))
                  AB=AB+AA
                END IF
              END IF
    6       CONTINUE
            IF(DABS(AB).GT.EPS) THEN
	      ABB=AB
              IF((ICOLOM+ISOTOP).EQ.1) THEN
                AB=A1*AB
                IF(DABS(AB).GT.EPS)
     :                 CALL SAVENON(2,AB,KL,0,LA,0,LB,JA,JB,0)
              END IF
C Orbit 1221
              IF(IORBORB.EQ.1) THEN
                IF(DABS(A2).GT.EPS) THEN
	          KLL=KL-1
                  ABB=A2*ABB*HALF/DSQRT(DBLE(2*KL+1))
                  CALL SAVENON(9,ABB,KLL,LA,LB,LB,LA,JA,JB,0)
                  CALL SAVENON(9,ABB,KLL,LB,LA,LA,LB,JA,JB,0)
C                 WRITE(79,656) ABB,KLL,LJ(IA),LJ(IB),JA,JB
  656 FORMAT(1X,'1221 O-O','AA=',F17.7,'K=',I3,'LA=',2I3,'JA=',I4,
     :'JB=',I4)
                END IF
              END IF
            END IF
          END IF
        END IF
    5 CONTINUE
      RETURN
      END
