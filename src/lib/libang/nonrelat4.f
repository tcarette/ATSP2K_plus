*
*     -------------------------------------------------------------
*      N O N R E L A T 4 
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :      N'1 = N1 +- 1        *
*                                             N'2 = N2 +- 1        *
*                                             N'3 = N3 -+ 2        *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*                                                                  *
      SUBROUTINE NONRELAT4(IA,IB,IC,ID)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON /OPERAT/ ICOLOM,ISOTOP,IORBORB
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      IF(IHSH.LE.2)RETURN
      IF(IA.EQ.IB) THEN
        IF(ISOTOP.EQ.1) THEN
	  IF(ITTK(LJ(IC),LJ(IA),1).EQ.0)RETURN
	  IF(ITTK(LJ(ID),LJ(IA),1).EQ.0)RETURN
        ENDIF
        CALL NONRELAT41(IC,ID,IA,1,IA,IB,IC,ID)
      ELSEIF(IC.EQ.ID) THEN
        IF(ISOTOP.EQ.1) THEN
	  IF(ITTK(LJ(IA),LJ(IC),1).EQ.0)RETURN
	  IF(ITTK(LJ(IB),LJ(IC),1).EQ.0)RETURN
        ENDIF
        CALL NONRELAT41(IA,IB,IC,2,IA,IB,IC,ID)
      ELSE
	WRITE(6,'(A)')  ' KLAIDA NONRELAT4  '
        STOP
      ENDIF
      RETURN
      END
