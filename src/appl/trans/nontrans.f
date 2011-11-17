*
*     --------------------------------------------------------------
*      N O N T R A N S
*     --------------------------------------------------------------
*
*     THE ROUTINE EVALUATES THE TRANSITION OPERATORS WITH          *
*     ORTHOGONAL ORBITALS                                          *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville           September 1997   * 
*
      SUBROUTINE NONTRANS(KA,KB,CL,CV)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL REL,VOK
      COMMON /EMS/IEM(4),IFL,JI,JF,LAM,REL,VOK
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/SAVECOM/LRHO,LSIG,CL2,CV2
      EXTERNAL TRANSITION
      CL2=ZERO
      CV2=ZERO
      CL=CL2
      CV=CV2
      IX=0
      IRHO=0
      IRHOP=0
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
      IF(IX.EQ.2) THEN
	LRHO=LJ(IRHO)
	LSIG=LJ(IRHOP)
	IF(IFL.EQ.1) THEN
	  IF(ITTK(LRHO,LSIG,LAM).EQ.0) RETURN
        ELSE
	  IF(ITTK(LRHO,LSIG,LAM-1).EQ.0) RETURN
	ENDIF
        CALL ONEPARTICLE2(KA,KB,IRHO,IRHOP,TRANSITION)
      ELSEIF(IX.EQ.0) THEN
        DO 2 K1=1,IHSH
          IF(NOSH1(K1).NE.0) THEN
	    LRHO=LJ(K1)
	    LSIG=LRHO
            IF(IFL.EQ.1) THEN
	      IF(ITTK(LRHO,LSIG,LAM).NE.0) 
     :             CALL ONEPARTICLE1(KA,KB,K1,TRANSITION)
            ELSE
	      IF(ITTK(LRHO,LSIG,LAM-1).NE.0) 
     :             CALL ONEPARTICLE1(KA,KB,K1,TRANSITION)
	    ENDIF
          ENDIF
    2   CONTINUE
      ENDIF
      CL=CL2
      CV=CV2
      RETURN
      END
