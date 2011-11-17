*
*     ------------------------------------------------------------------
*	C F G N 1
*     ------------------------------------------------------------------
*
*	Read the configurations for a state and determine the
*	non-orthogonal orbitals
*
      SUBROUTINE CFGN1(INPUT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NCD=2000)
*
      CHARACTER BUFFER*3
      CHARACTER*1 JAJCMP(NWD,3), INPUT*24
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,JSC0,JSC1,
     :JSC2,JSC3
      COMMON/OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
     : ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP,
     :     IORTH(NWD*(NWD-1)/2)
      COMMON/STATES/NCFG,MAXORB,IAJCMP(NWD),LJCOMP(NWD),
     :NJCOMP(NWD),NOCCSH(NCD),NELCSH(5,NCD),NOCORB(5,NCD),J1QNRD(9,NCD)
      INTEGER IEL(2),IBUFF(2)
      DATA IASTER,IBLANK/3H*  ,3H   /
*
      CALL CFGO1(NCFG,MAXORB,IAJCMP,LJCOMP,NJCOMP,NOCCSH,NELCSH,
     :            NOCORB,J1QNRD,NCD,INPUT)
*
*  ---  SEPARATE THE ELECTRON LABEL CHARACTERS AND LEFT JUSTIFY
*
      DO 10 I = 1,MAXORB
         WRITE(BUFFER,'(A3)') IAJCMP(I)
	 READ(BUFFER,'(3A1)') (JAJCMP(I,J),J=1,3)
	 IF (JAJCMP(I,1) .EQ. ' ') THEN
	    JAJCMP(I,1) = JAJCMP(I,2)
	    JAJCMP(I,2) = JAJCMP(I,3)
	    JAJCMP(I,3) = ' '
	 END IF
10    CONTINUE
*
*  ---  INITIALIZE THE ORTHOGONALITY ARRAY
*
      M1 = (MAXORB*(MAXORB-1))/2
      DO 20 I = 1,M1
         IORTH(I) = -1
20    CONTINUE
*
*  ---  SET ORBITALS IN THE SAME CONFIGURATION TO BE ORTHOGONAL
*
      DO 63 I = 1,NCFG
      N = NOCCSH(I)
      DO 66 J = 1,N-1
      DO 66 JJ = J+1,N
         I1 = NOCORB(J,I)
         J1 = NOCORB(JJ,I)
         IF (J1 .GT. I1) THEN
            M = I1
            I1 = J1
            J1 = M
         ENDIF
         IORTH(J1+((I1-1)*(I1-2))/2) = 0
   66 CONTINUE
   63 CONTINUE
*
* --- DETERMINE THE NON-ORTHOGONAL ORBITALS
*
      NORTH = 0
      DO 18 J = 1,MAXORB-1
      DO 19 I = J+1,MAXORB
         IJ = J + ((I-1)*(I-2))/2
         IF (JAJCMP(I,2) .EQ. JAJCMP(J,2) .AND.
     :       JAJCMP(I,3) .NE. ' ' .AND.
     :       JAJCMP(J,3) .NE. ' ' .AND.
     :       JAJCMP(I,3) .NE. JAJCMP(J,3) .AND.
     :       IORTH(IJ) .NE. 0 ) THEN
             NORTH = NORTH + 1
             IORTH(IJ) = 1
         ENDIF
19    CONTINUE
18    CONTINUE
*
      READ(IREAD,*,END=90)
 79   READ(IREAD,'(2(1X,A3))',END=90) IBUFF(1),IBUFF(2)
      IF (IBUFF(1) .NE. IASTER .AND. IBUFF(1) .NE. IBLANK) THEN
         DO 80 I = 1,2
            DO 81 J = 1,MAXORB
               IF (IBUFF(I) .EQ. IAJCMP(J)) THEN
                  IEL(I) = J
                  GO TO 80
               ENDIF
 81         CONTINUE
            WRITE(*,'(A,A3,A)') ' ELECTRON ',IBUFF(I),' NOT FOUND'
            STOP
 80      CONTINUE
         IF (IEL(1) .GT. IEL(2) ) THEN
            I = IEL(1)
            IEL(1) = IEL(2)
            IEL(2) = I
         END IF
         IJ = IEL(1) + ((IEL(2)-1)*(IEL(2)-2))/2
         IF (IORTH(IJ) .EQ. 1) NORTH = NORTH - 1
         IORTH(IJ) = 0
         WRITE(IWRITE,'(1X,A3,A,A3)')
     :           IBUFF(1),' is orthogonal to ',IBUFF(2)
         GO TO 79
      END IF
  90  RETURN
      END
