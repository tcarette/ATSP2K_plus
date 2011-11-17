*
*     --------------------------------------------------------------
*      N O N H I P E R
*     --------------------------------------------------------------
*
*     THE ROUTINE EVALUATES THE TRANSITION OPERATORS WITH          *
*     ORTHOGONAL ORBITALS                                          *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville           September 1997   * 
*
      SUBROUTINE NONHIPER(IREZ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/HYPER/VHY(20),NHY,IRHY(20),ISHY(20),JRHY(20),JSHY(20),
     1 IMUHY(20),IMUPHY(20),JMUHY(20),JMUPHY(20),
     2 INUHY(20),INUPHY(20),JNUHY(20),JNUPHY(20)
      COMMON/PAPIL/ IIRHO,IISIG
      EXTERNAL ORBIT,DIPOL,CONTAC,QDRPOL,VOLUM
      NHY=0
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
	IIRHO=IRHO
	IISIG=IRHOP
        IF(IREZ.EQ.1) THEN
*         ORBITA
          CALL ONEPARTICLE2(1,0,IRHO,IRHOP,ORBIT)
        ELSEIF(IREZ.EQ.2) THEN
*         DIPOL
          CALL ONEPARTICLE2(2,1,IRHO,IRHOP,DIPOL)
        ELSEIF(IREZ.EQ.3) THEN
*         CONTACT
          CALL ONEPARTICLE2(0,1,IRHO,IRHOP,CONTAC)
        ELSEIF(IREZ.EQ.4) THEN
*         QDRPOLE
          CALL ONEPARTICLE2(2,0,IRHO,IRHOP,QDRPOL)
        ELSEIF(IREZ.EQ.5) THEN
*         VOLUM
          CALL ONEPARTICLE2(0,0,IRHO,IRHOP,VOLUM)
        ELSE
          RETURN
        ENDIF

      ELSEIF(IX.EQ.0) THEN
        DO 2 K1=1,IHSH
          IF(NOSH1(K1).NE.0) THEN
	    LRHO=LJ(K1)
	    LSIG=LRHO
	    IIRHO=K1
	    IISIG=K1
            IF(IREZ.EQ.1) THEN
*             ORBITA
              CALL ONEPARTICLE1(1,0,K1,ORBIT)
            ELSEIF(IREZ.EQ.2) THEN
*             DIPOL
              CALL ONEPARTICLE1(2,1,K1,DIPOL)
            ELSEIF(IREZ.EQ.3) THEN
*             CONTACT
              CALL ONEPARTICLE1(0,1,K1,CONTAC)
            ELSEIF(IREZ.EQ.4) THEN
*             QDRPOLE
              CALL ONEPARTICLE1(2,0,K1,QDRPOL)
            ELSEIF(IREZ.EQ.5) THEN
*             VOLUM
              CALL ONEPARTICLE1(0,0,K1,VOLUM)
            ELSE
              RETURN
            ENDIF
          ENDIF
    2   CONTINUE
      ENDIF
      RETURN
      END
