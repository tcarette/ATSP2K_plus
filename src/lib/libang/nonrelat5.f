*
*     -------------------------------------------------------------
*      N O N R E L A T 5 
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :    N'1 = N1 (+-) 1        *
*                                           N'2 = N2 (+-) 1        *
*                                           N'3 = N3 (+-) 1        *
*                                           N'4 = N4 (+-) 1        *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*
      SUBROUTINE NONRELAT5(IA,IB,IC,ID)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      IF(IHSH.LE.3)RETURN
      IF(IB.LT.IC) THEN
        CALL NONRELAT51(IA,IB,IC,ID,1)
      ELSEIF(IA.GT.ID.AND.IB.GT.ID) THEN
        CALL NONRELAT51(IC,ID,IA,IB,2)
      ELSEIF(IB.GT.IC.AND.IB.LT.ID.AND.IA.LT.IC) THEN
        CALL NONRELAT52(IA,IC,IB,ID,1)
      ELSEIF(IB.GT.IC.AND.IB.GT.ID.AND.IA.GT.IC) THEN
        CALL NONRELAT52(IC,IA,ID,IB,2)
      ELSEIF(IB.GT.IC.AND.IB.GT.ID.AND.IA.LT.IC) THEN
        CALL NONRELAT53(IA,IC,ID,IB,1)
      ELSEIF(IB.GT.IC.AND.IB.LT.ID.AND.IA.GT.IC) THEN
        CALL NONRELAT53(IC,IA,IB,ID,2)
      ELSE
	WRITE(6,'(A)')  ' KLAIDA NONRELAT5  '
      STOP
      ENDIF
      RETURN
      END
