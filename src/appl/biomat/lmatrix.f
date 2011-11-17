*
*     ------------------------------------------------------------------
*      L M A T R I X
*     ------------------------------------------------------------------
*
*     THE ROUTINE EVALUATES THE ONE-ELECTRON NON-RELATIVISTIC 
*     HAMILTONIAN WITH ORTHOGONAL ORBITALS
*
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville           September 1997   * 
*
      SUBROUTINE LMATRIX
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/DIAGNL/IDIAG,JA,JB
      IX=0
      DO 1 J=1,IHSH
        N=NOSH1(J)-NOSH2(J)
        IF(IABS(N).GT.1) RETURN
        IF(N.EQ.1) THEN
          IRHO = J
          IX =IX+1
        ELSEIF(N+1.EQ.0) THEN
          IRHOP = J
          IX=IX+1
        ENDIF
    1 CONTINUE
      IF(IX.GT.2) RETURN
      N1=2*IHSH-1
      IF((J1QN1(N1,2)-J1QN2(N1,2)).NE.0)RETURN
      IF((J1QN1(N1,3)-J1QN2(N1,3)).NE.0)RETURN
      IF(IX.EQ.2) THEN
	IF(LJ(IRHOP).NE.LJ(IRHO)) RETURN
        CALL LMATRIX2(IRHOP,IRHO)
      ELSEIF(IX.EQ.0) THEN
        DO 2 K1=1,IHSH
          IF(NOSH1(K1).NE.0) CALL LMATRIX1(K1)
    2   CONTINUE
      ENDIF
      RETURN
      END
