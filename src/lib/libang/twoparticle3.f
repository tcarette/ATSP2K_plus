*
*     -------------------------------------------------------------
*      T W O P A R T I C L E 3 
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
*                                              N'2 = N2 + 1        *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE TWOPARTICLE3(IA,IB,IC,ID)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON /CASEOP/ IOCASE
      IF(IHSH.LE.1)RETURN
*
      IF(IOCASE.EQ.1) THEN
* Spin - spin      operator
	KGAL=2
      ELSEIF(IOCASE.EQ.2) THEN
* Spin-other-orbit operator
	KGAL=1
      ENDIF
*
      IF(IB.EQ.ID) THEN
        IF(IA.EQ.IB.OR.IC.EQ.IB) THEN
          IF(IA.EQ.IC)GO TO 10
          IIA=MIN0(IA,IC)
          IIB=MAX0(IA,IC)
          CALL RLSP0(1,IIA,IIB,0,IAT)
          IF(IAT.EQ.0) RETURN
          CALL RLSP0(2,IIA,IIB,2*KGAL,IAT)
          IF(IAT.EQ.0) RETURN
          CALL RLSP0(3,IIA,IIB,2*KGAL,IAT)
          IF(IAT.EQ.0) RETURN
          IF(IC.EQ.IB) THEN
            CALL TWOPARTICLE31(IC,IA,IA,IB,IC,ID)
          ELSE
            CALL TWOPARTICLE32(IC,IA,IA,IB,IC,ID)
	  ENDIF
	ELSE
          CALL EILE(IA,IB,IC,IAA,IBB,ICC)
          CALL RLSP00(1,IAA,IBB,ICC,ICC,0,IAT)
          IF(IAT.EQ.0) RETURN
          CALL RLSP00(2,IAA,IBB,ICC,ICC,2*KGAL,IAT)
          IF(IAT.EQ.0) RETURN
          CALL RLSP00(3,IAA,IBB,ICC,ICC,2*KGAL,IAT)
          IF(IAT.EQ.0) RETURN
          CALL TWOPARTICLE33(IC,IA,IB,IA,IB,IC,ID)
	ENDIF
      ELSEIF(IA.EQ.IC) THEN
        IF(IB.EQ.IA.OR.ID.EQ.IA) THEN
          IF(IB.EQ.ID)GO TO 10
          IIA=MIN0(IB,ID)
          IIB=MAX0(IB,ID)
          CALL RLSP0(1,IIA,IIB,0,IAT)
          IF(IAT.EQ.0) RETURN
          CALL RLSP0(2,IIA,IIB,2*KGAL,IAT)
          IF(IAT.EQ.0) RETURN
          CALL RLSP0(3,IIA,IIB,2*KGAL,IAT)
          IF(IAT.EQ.0) RETURN
          IF(ID.EQ.IA) THEN
            CALL TWOPARTICLE31(ID,IB,IA,IB,IC,ID)
          ELSE
            CALL TWOPARTICLE32(ID,IB,IA,IB,IC,ID)
	  ENDIF
	ELSE
          CALL EILE(IA,IB,ID,IAA,IBB,ICC)
          CALL RLSP00(1,IAA,IBB,ICC,ICC,0,IAT)
          IF(IAT.EQ.0) RETURN
          CALL RLSP00(2,IAA,IBB,ICC,ICC,2*KGAL,IAT)
          IF(IAT.EQ.0) RETURN
          CALL RLSP00(3,IAA,IBB,ICC,ICC,2*KGAL,IAT)
          IF(IAT.EQ.0) RETURN
          CALL TWOPARTICLE33(ID,IB,IA,IA,IB,IC,ID)
	ENDIF
      ELSEIF(IA.EQ.ID) THEN
        IF(IB.EQ.IA.OR.IC.EQ.IA) THEN
          IF(IB.EQ.IC)GO TO 10
          IIA=MIN0(IB,IC)
          IIB=MAX0(IB,IC)
          CALL RLSP0(1,IIA,IIB,0,IAT)
          IF(IAT.EQ.0) RETURN
          CALL RLSP0(2,IIA,IIB,2*KGAL,IAT)
          IF(IAT.EQ.0) RETURN
          CALL RLSP0(3,IIA,IIB,2*KGAL,IAT)
          IF(IAT.EQ.0) RETURN
          IF(IC.EQ.ID) THEN
            CALL TWOPARTICLE31(IC,IB,IA,IB,IC,ID)
          ELSE
            CALL TWOPARTICLE32(IC,IB,IA,IB,IC,ID)
	  ENDIF
	ELSE
          CALL EILE(IA,IB,IC,IAA,IBB,ICC)
          CALL RLSP00(1,IAA,IBB,ICC,ICC,0,IAT)
          IF(IAT.EQ.0) RETURN
          CALL RLSP00(2,IAA,IBB,ICC,ICC,2*KGAL,IAT)
          IF(IAT.EQ.0) RETURN
          CALL RLSP00(3,IAA,IBB,ICC,ICC,2*KGAL,IAT)
          IF(IAT.EQ.0) RETURN
          CALL TWOPARTICLE33(IC,IB,IA,IA,IB,ID,IC)
	ENDIF
      ELSEIF(IB.EQ.IC) THEN
        IF(IA.EQ.IB.OR.ID.EQ.IB) THEN
          IF(IA.EQ.ID)GO TO 10
          IIA=MIN0(IA,ID)
          IIB=MAX0(IA,ID)
          CALL RLSP0(1,IIA,IIB,0,IAT)
          IF(IAT.EQ.0) RETURN
          CALL RLSP0(2,IIA,IIB,2*KGAL,IAT)
          IF(IAT.EQ.0) RETURN
          CALL RLSP0(3,IIA,IIB,2*KGAL,IAT)
          IF(IAT.EQ.0) RETURN
          IF(ID.EQ.IB) THEN
            CALL TWOPARTICLE31(ID,IA,IA,IB,IC,ID)
          ELSE
            CALL TWOPARTICLE32(ID,IA,IA,IB,IC,ID)
	  ENDIF
	ELSE
          CALL EILE(IA,IB,ID,IAA,IBB,ICC)
          CALL RLSP00(1,IAA,IBB,ICC,ICC,0,IAT)
          IF(IAT.EQ.0) RETURN
          CALL RLSP00(2,IAA,IBB,ICC,ICC,2*KGAL,IAT)
          IF(IAT.EQ.0) RETURN
          CALL RLSP00(3,IAA,IBB,ICC,ICC,2*KGAL,IAT)
          IF(IAT.EQ.0) RETURN
          CALL TWOPARTICLE33(ID,IA,IB,IA,IB,ID,IC)
	ENDIF
      ENDIF
      RETURN
   10 WRITE(6,'(A)') ' ERROR IN NONRELAT3'
      STOP
      END
