*
*     -------------------------------------------------------------
*      T W O P A R T I C L E 5
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :    N'1 = N1 (+-) 1        *
*                                           N'2 = N2 (+-) 1        *
*                                           N'3 = N3 (+-) 1        *
*                                           N'4 = N4 (+-) 1        *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE TWOPARTICLE5(IA,IB,IC,ID)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      COMMON/TRK2/BD3(3),BD4(3),BK3(3),BK4(3),
     *ID3(7),ID4(7),IK3(7),IK4(7)
      COMMON/PERMAT/IRS(2,2),IRL(20,20),RS(2,2),RL(20,20)
      COMMON/KAMPAS/ IW1(2,20),IW2(2,20),IWAA(2,2,20,20),
     :RW1(2,20),RW2(2,20),RWAA(2,2,20,20)
      COMMON /CASEOP/ IOCASE
      EXTERNAL TWO51,TWO52,TWO53,TWO54,TWO55,TWO56
      IF(IHSH.LE.3)RETURN
*
      IF(IOCASE.EQ.1) THEN
* Spin - spin      operator
	KGAL=2
      ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
	KGAL=1
      ENDIF
      DO 13 I1=1,2
        DO 14 I2=1,2
          IRS(I1,I2)=0
   14   CONTINUE
   13 CONTINUE
      DO 15 I1=1,20
        IW1(1,I1)=0
        IW1(2,I1)=0
        IW2(1,I1)=0
        IW2(2,I1)=0
        DO 16 I2=1,20
          IRL(I1,I2)=0
          IWAA(1,1,I1,I2)=0
          IWAA(2,1,I1,I2)=0
   16   CONTINUE
   15 CONTINUE
*   i)
      IF(IB.LT.IC) THEN
        CALL RLSP00(1,IA,IB,IC,ID,0,IAT)
        IF(IAT.EQ.0) RETURN
        CALL RLSP00(2,IA,IB,IC,ID,2*KGAL,IAT)
        IF(IAT.EQ.0) RETURN
        CALL RLSP00(3,IA,IB,IC,ID,2*KGAL,IAT)
        IF(IAT.EQ.0) RETURN
        CALL A1A2A3A4LSP(IA,IB,IC,ID,HALF,HALF,-HALF,-HALF,W)
        IF(DABS(W).LT.EPS) RETURN 
        IF(IOCASE.EQ.1) THEN
* Spin-Spin operator
          IP=ITREXG2(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK)+1
        ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
          IP=ITREXG(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK)+1
        ENDIF
        IF(IKK.GT.0) THEN
          IG=IP+IKK-1
          DO 1 I=IP,IG
            KL=I-1
*
            IF(IOCASE.EQ.1) THEN
* Spin-spin operator
              CALL SSC(IG,KL,IA,IB,IC,ID,IA,IB,IC,ID,1,TWO51)
            ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
              CALL SOOC(IG,KL,IA,IB,IC,ID,IA,IB,IC,ID,1,TWO51)
            ENDIF
*
    1     CONTINUE
        ENDIF
        IF(IOCASE.EQ.1) THEN
* Spin-Spin operator
          IP=ITREXG2(LJ(IA),LJ(ID),LJ(IB),LJ(IC),IKK)+1
        ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
          IP=ITREXG(LJ(IA),LJ(ID),LJ(IB),LJ(IC),IKK)+1
        ENDIF
        IF(IKK.LE.0)RETURN
        IG=IP+IKK-1
        DO 2 I=IP,IG
          KL=I-1
*
          IF(IOCASE.EQ.1) THEN
* Spin-spin operator
            CALL SSC(IG,KL,IA,IB,IC,ID,IA,IB,ID,IC,1,TWO52)
          ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
            CALL SOOC(IG,KL,IA,IB,IC,ID,IA,IB,ID,IC,1,TWO52)
          ENDIF
*
    2   CONTINUE
*
      ELSEIF(IA.GT.ID.AND.IB.GT.ID) THEN
        CALL RLSP00(1,IC,ID,IA,IB,0,IAT)
        IF(IAT.EQ.0) RETURN
        CALL RLSP00(2,IC,ID,IA,IB,2*KGAL,IAT)
        IF(IAT.EQ.0) RETURN
        CALL RLSP00(3,IC,ID,IA,IB,2*KGAL,IAT)
        IF(IAT.EQ.0) RETURN
        CALL A1A2A3A4LSP(IC,ID,IA,IB,-HALF,-HALF,HALF,HALF,W)
        IF(DABS(W).LT.EPS) RETURN 
        IF(IOCASE.EQ.1) THEN
* Spin-Spin operator
          IP=ITREXG2(LJ(IC),LJ(IA),LJ(ID),LJ(IB),IKK)+1
        ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
          IP=ITREXG(LJ(IC),LJ(IA),LJ(ID),LJ(IB),IKK)+1
        ENDIF
        IF(IKK.GT.0) THEN
          IG=IP+IKK-1
          DO 3 I=IP,IG
            KL=I-1
*
            IF(IOCASE.EQ.1) THEN
* Spin-spin operator
              CALL SSC(IG,KL,IC,ID,IA,IB,IA,IB,IC,ID,2,TWO51)
            ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
              CALL SOOC(IG,KL,IC,ID,IA,IB,IA,IB,IC,ID,2,TWO51)
            ENDIF
*
    3     CONTINUE
        ENDIF
        IF(IOCASE.EQ.1) THEN
* Spin-Spin operator
          IP=ITREXG2(LJ(IC),LJ(IB),LJ(ID),LJ(IA),IKK)+1
        ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
          IP=ITREXG(LJ(IC),LJ(IB),LJ(ID),LJ(IA),IKK)+1
        ENDIF
        IF(IKK.LE.0)RETURN
        IG=IP+IKK-1
        DO 4 I=IP,IG
          KL=I-1
*
          IF(IOCASE.EQ.1) THEN
* Spin-spin operator
            CALL SSC(IG,KL,IC,ID,IA,IB,IB,IA,IC,ID,2,TWO52)
          ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
            CALL SOOC(IG,KL,IC,ID,IA,IB,IB,IA,IC,ID,2,TWO52)
          ENDIF
*
    4   CONTINUE
*  ii)
      ELSEIF(IB.GT.IC.AND.IB.LT.ID.AND.IA.LT.IC) THEN
        CALL RLSP00(1,IA,IC,IB,ID,0,IAT)
        IF(IAT.EQ.0) RETURN
        CALL RLSP00(2,IA,IC,IB,ID,2*KGAL,IAT)
        IF(IAT.EQ.0) RETURN
        CALL RLSP00(3,IA,IC,IB,ID,2*KGAL,IAT)
        IF(IAT.EQ.0) RETURN
        CALL A1A2A3A4LSP(IA,IC,IB,ID,HALF,-HALF,HALF,-HALF,W)
        IF(DABS(W).LT.EPS) RETURN 
        IF(IOCASE.EQ.1) THEN
* Spin-Spin operator
          IP=ITREXG2(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK)+1
        ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
          IP=ITREXG(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK)+1
        ENDIF
        IF(IKK.GT.0) THEN
          IG=IP+IKK-1
          DO 5 I=IP,IG
            KL=I-1
*
            IF(IOCASE.EQ.1) THEN
* Spin-spin operator
              CALL SSC(IG,KL,IA,IC,IB,ID,IA,IB,IC,ID,1,TWO53)
            ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
              CALL SOOC(IG,KL,IA,IC,IB,ID,IA,IB,IC,ID,1,TWO53)
            ENDIF
*
    5     CONTINUE
        ENDIF
        IF(IOCASE.EQ.1) THEN
* Spin-Spin operator
          IP=ITREXG2(LJ(IA),LJ(ID),LJ(IB),LJ(IC),IKK)+1
        ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
          IP=ITREXG(LJ(IA),LJ(ID),LJ(IB),LJ(IC),IKK)+1
        ENDIF
        IF(IKK.LE.0)RETURN
        IG=IP+IKK-1
        DO 6 I=IP,IG
          KL=I-1
*
          IF(IOCASE.EQ.1) THEN
* Spin-spin operator
            CALL SSC(IG,KL,IA,IC,IB,ID,IA,IB,ID,IC,1,TWO54)
          ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
            CALL SOOC(IG,KL,IA,IC,IB,ID,IA,IB,ID,IC,1,TWO54)
          ENDIF
*
    6   CONTINUE
*
      ELSEIF(IB.GT.IC.AND.IB.GT.ID.AND.IA.GT.IC) THEN
        CALL RLSP00(1,IC,IA,ID,IB,0,IAT)
        IF(IAT.EQ.0) RETURN
        CALL RLSP00(2,IC,IA,ID,IB,2*KGAL,IAT)
        IF(IAT.EQ.0) RETURN
        CALL RLSP00(3,IC,IA,ID,IB,2*KGAL,IAT)
        IF(IAT.EQ.0) RETURN
        CALL A1A2A3A4LSP(IC,IA,ID,IB,-HALF,HALF,-HALF,HALF,W)
        IF(DABS(W).LT.EPS) RETURN 
        IF(IOCASE.EQ.1) THEN
* Spin-Spin operator
          IP=ITREXG2(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK)+1
        ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
          IP=ITREXG(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK)+1
        ENDIF
        IF(IKK.GT.0) THEN
          IG=IP+IKK-1
          DO 7 I=IP,IG
            KL=I-1
*
            IF(IOCASE.EQ.1) THEN
* Spin-spin operator
              CALL SSC(IG,KL,IC,IA,ID,IB,IA,IB,IC,ID,2,TWO53)
            ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
              CALL SOOC(IG,KL,IC,IA,ID,IB,IA,IB,IC,ID,2,TWO53)
            ENDIF
*
    7     CONTINUE
        ENDIF
        IF(IOCASE.EQ.1) THEN
* Spin-Spin operator
          IP=ITREXG2(LJ(IA),LJ(ID),LJ(IB),LJ(IC),IKK)+1
        ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
          IP=ITREXG(LJ(IA),LJ(ID),LJ(IB),LJ(IC),IKK)+1
        ENDIF
        IF(IKK.LE.0)RETURN
        IG=IP+IKK-1
        DO 8 I=IP,IG
          KL=I-1
*
          IF(IOCASE.EQ.1) THEN
* Spin-spin operator
            CALL SSC(IG,KL,IC,IA,ID,IB,IB,IA,IC,ID,2,TWO54)
          ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
            CALL SOOC(IG,KL,IC,IA,ID,IB,IB,IA,IC,ID,2,TWO54)
          ENDIF
*
    8   CONTINUE
* iii)
      ELSEIF(IB.GT.IC.AND.IB.GT.ID.AND.IA.LT.IC) THEN
        CALL RLSP00(1,IA,IC,ID,IB,0,IAT)
        IF(IAT.EQ.0) RETURN
        CALL RLSP00(2,IA,IC,ID,IB,2*KGAL,IAT)
        IF(IAT.EQ.0) RETURN
        CALL RLSP00(3,IA,IC,ID,IB,2*KGAL,IAT)
        IF(IAT.EQ.0) RETURN
        CALL A1A2A3A4LSP(IA,IC,ID,IB,HALF,-HALF,-HALF,HALF,W)
        IF(DABS(W).LT.EPS) RETURN 
        IF(IOCASE.EQ.1) THEN
* Spin-Spin operator
          IP=ITREXG2(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK)+1
        ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
          IP=ITREXG(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK)+1
        ENDIF
        IF(IKK.GT.0) THEN
          IG=IP+IKK-1
          DO 9 I=IP,IG
            KL=I-1
*
            IF(IOCASE.EQ.1) THEN
* Spin-spin operator
              CALL SSC(IG,KL,IA,IC,ID,IB,IA,IB,IC,ID,1,TWO55)
            ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
              CALL SOOC(IG,KL,IA,IC,ID,IB,IA,IB,IC,ID,1,TWO55)
            ENDIF
*
    9     CONTINUE
        ENDIF
        IF(IOCASE.EQ.1) THEN
* Spin-Spin operator
          IP=ITREXG2(LJ(IA),LJ(ID),LJ(IB),LJ(IC),IKK)+1
        ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
          IP=ITREXG(LJ(IA),LJ(ID),LJ(IB),LJ(IC),IKK)+1
        ENDIF
        IF(IKK.LE.0)RETURN
        IG=IP+IKK-1
        DO 10 I=IP,IG
          KL=I-1
*
          IF(IOCASE.EQ.1) THEN
* Spin-spin operator
            CALL SSC(IG,KL,IA,IC,ID,IB,IA,IB,ID,IC,1,TWO56)
          ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
            CALL SOOC(IG,KL,IA,IC,ID,IB,IA,IB,ID,IC,1,TWO56)
          ENDIF
*
   10   CONTINUE
*
      ELSEIF(IB.GT.IC.AND.IB.LT.ID.AND.IA.GT.IC) THEN
        CALL RLSP00(1,IC,IA,IB,ID,0,IAT)
        IF(IAT.EQ.0) RETURN
        CALL RLSP00(2,IC,IA,IB,ID,2*KGAL,IAT)
        IF(IAT.EQ.0) RETURN
        CALL RLSP00(3,IC,IA,IB,ID,2*KGAL,IAT)
        IF(IAT.EQ.0) RETURN
        CALL A1A2A3A4LSP(IC,IA,IB,ID,-HALF,HALF,HALF,-HALF,W)
        IF(DABS(W).LT.EPS) RETURN 
        IF(IOCASE.EQ.1) THEN
* Spin-Spin operator
          IP=ITREXG2(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK)+1
        ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
          IP=ITREXG(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK)+1
        ENDIF
        IF(IKK.GT.0) THEN
          IG=IP+IKK-1
          DO 11 I=IP,IG
            KL=I-1
*
            IF(IOCASE.EQ.1) THEN
* Spin-spin operator
              CALL SSC(IG,KL,IC,IA,IB,ID,IA,IB,IC,ID,2,TWO55)
            ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
              CALL SOOC(IG,KL,IC,IA,IB,ID,IA,IB,IC,ID,2,TWO55)
            ENDIF
*
   11     CONTINUE
        ENDIF
        IF(IOCASE.EQ.1) THEN
* Spin-Spin operator
          IP=ITREXG2(LJ(IA),LJ(ID),LJ(IB),LJ(IC),IKK)+1
        ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
          IP=ITREXG(LJ(IA),LJ(ID),LJ(IB),LJ(IC),IKK)+1
        ENDIF
        IF(IKK.LE.0)RETURN
        IG=IP+IKK-1
        DO 12 I=IP,IG
          KL=I-1
*
          IF(IOCASE.EQ.1) THEN
* Spin-spin operator
            CALL SSC(IG,KL,IC,IA,IB,ID,IB,IA,IC,ID,2,TWO56)
          ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
            CALL SOOC(IG,KL,IC,IA,IB,ID,IB,IA,IC,ID,2,TWO56)
          ENDIF
*
   12   CONTINUE
*
      ELSE
	WRITE(6,'(A)')  ' ERROR IN SUBROUTINE  TWOPARTICLE5  '
      STOP
      ENDIF
      RETURN
      END
