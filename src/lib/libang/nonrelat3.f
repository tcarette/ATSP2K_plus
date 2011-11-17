*
*     -------------------------------------------------------------
*      N O N R E L A T 3 
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
*                                              N'2 = N2 + 1        *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*
      SUBROUTINE NONRELAT3(IA,IB,IC,ID,IIRE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON /OPERAT/ ICOLOM,ISOTOP,IORBORB
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      IF(IHSH.LE.1)RETURN
      IF(IB.EQ.ID) THEN
        IF(IA.EQ.IB.OR.IC.EQ.IB) THEN
          IF(IA.EQ.IC)GO TO 10
          IF(IC.EQ.IB) THEN
            CALL NONRELAT31(IC,IA,IA,IB,IC,ID,IIRE)
          ELSE
            IF(IIRE.EQ.0)RETURN
            IF(ISOTOP.EQ.1) THEN
	      IF(ITTK(LJ(IA),LJ(IA),1).EQ.0)RETURN
	      IF(ITTK(LJ(IC),LJ(IA),1).EQ.0)RETURN
            ENDIF
            CALL NONRELAT32(IC,IA,IA,IB,IC,ID)
	  ENDIF
	ELSE
          IF(IIRE.EQ.0)RETURN
          CALL NONRELAT33(IC,IA,IB,1,IA,IB,IC,ID)
	ENDIF
      ELSEIF(IA.EQ.IC) THEN
        IF(IB.EQ.IA.OR.ID.EQ.IA) THEN
          IF(IB.EQ.ID)GO TO 10
          IF(ID.EQ.IA) THEN
           CALL NONRELAT31(ID,IB,IA,IB,IC,ID,IIRE)
          ELSE
            IF(IIRE.EQ.0)RETURN
            IF(ISOTOP.EQ.1) THEN
	      IF(ITTK(LJ(IB),LJ(IB),1).EQ.0)RETURN
	      IF(ITTK(LJ(ID),LJ(IB),1).EQ.0)RETURN
            ENDIF
            CALL NONRELAT32(ID,IB,IA,IB,IC,ID)
	  ENDIF
	ELSE
          IF(IIRE.EQ.0)RETURN
          CALL NONRELAT33(ID,IB,IA,1,IA,IB,IC,ID)
	ENDIF
      ELSEIF(IA.EQ.ID) THEN
        IF(IB.EQ.IA.OR.IC.EQ.IA) THEN
          IF(IB.EQ.IC)GO TO 10
          IF(IC.EQ.ID) THEN
            CALL NONRELAT31(IC,IB,IA,IB,IC,ID,IIRE)
          ELSE
            IF(IIRE.EQ.0)RETURN
            IF(ISOTOP.EQ.1) THEN
	      IF(ITTK(LJ(IB),LJ(IB),1).EQ.0)RETURN
	      IF(ITTK(LJ(IC),LJ(IB),1).EQ.0)RETURN
            ENDIF
            CALL NONRELAT32(IC,IB,IA,IB,IC,ID)
	  ENDIF
	ELSE
          IF(IIRE.EQ.0)RETURN
          CALL NONRELAT33(IC,IB,IA,2,IA,IB,ID,IC)
	ENDIF
      ELSEIF(IB.EQ.IC) THEN
        IF(IA.EQ.IB.OR.ID.EQ.IB) THEN
          IF(IA.EQ.ID)GO TO 10
          IF(ID.EQ.IB) THEN
            CALL NONRELAT31(ID,IA,IA,IB,IC,ID,IIRE)
          ELSE
            IF(IIRE.EQ.0)RETURN
            IF(ISOTOP.EQ.1) THEN
	      IF(ITTK(LJ(IA),LJ(IA),1).EQ.0)RETURN
	      IF(ITTK(LJ(ID),LJ(IA),1).EQ.0)RETURN
            ENDIF
            CALL NONRELAT32(ID,IA,IA,IB,IC,ID)
	  ENDIF
	ELSE
          IF(IIRE.EQ.0)RETURN
          CALL NONRELAT33(ID,IA,IB,2,IA,IB,ID,IC)
	ENDIF
      ENDIF
      RETURN
   10 WRITE(6,'(A)') ' ERRO IN NONRELAT3'
      STOP
      END
