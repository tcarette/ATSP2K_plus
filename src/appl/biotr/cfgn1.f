*
*     ------------------------------------------------------------------
*       C F G N 1
*     ------------------------------------------------------------------
*
      SUBROUTINE CFGN1(NCLOSD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NWD=128,NWCD=20)
*
*       Read the configurations for a state and determine the
*       non-orthogonal orbitals
*
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,
     : JSC0,JSC(3),ISCW
      INTEGER IEL(2),IBUFF(2)
      CHARACTER JAJCMP(3,NWD)*1, BUFFER*8
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      POINTER(QIORTH,IORTH(1))
      COMMON/OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
     : ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP,
     :     QIORTH
      DATA iaster,iblnk/4h*   , 4h    /
*
      CALL CFGO1 (NCFG,MAXORB,NCLOSD,QIAJCMP,QLJCOMP,QNJCOMP,
     :                 QNOC,QNELCSH,QNOCORB,QJ1,QIAJCLD,QLJCLSD)
      NWF = MAXORB
 
*
*     Now deal with the non-orthogonality
*
      if (NWF .GT. 1) then
        print*,' alloc, qiorth: nwf*(nwf-1)/2 = ',nwf*(nwf-1)/2
        call alloc(qiorth,NWF*(NWF-1)/2,8)
      end if
*
*  ---  SEPARATE THE ELECTRON LABEL CHARACTERS AND LEFT JUSTIFY
*
       DO 10 I = 1,MAXORB
          WRITE(BUFFER,'(A3)') IAJCMP(I)
          READ(BUFFER,'(3A1)') (JAJCMP(J,I),J=1,3)
          IF (JAJCMP(1,I) .EQ. ' ') THEN
             JAJCMP(1,I) = JAJCMP(2,I)
             JAJCMP(2,I) = JAJCMP(3,I)
             JAJCMP(3,I) = ' '
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
      print*,' NCFG = ', NCFG
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
         IF (JAJCMP(2,I) .EQ. JAJCMP(2,J) .AND.
     :       JAJCMP(3,I) .NE. ' ' .AND.
     :       JAJCMP(3,J) .NE. ' ' .AND.
     :       JAJCMP(3,I) .NE. JAJCMP(3,J) .AND.
     :       IORTH(IJ) .NE. 0 ) THEN
             NORTH = NORTH + 1
             IORTH(IJ) = 1
         ENDIF
19    CONTINUE
18    CONTINUE
*
      READ(IREAD,*,END=90)
 79   READ(IREAD,'(2(1X,A3))',END=90) IBUFF(1),IBUFF(2)
      IF (IBUFF(1) .NE. iaster .AND. IBUFF(1) .NE. iblnk ) THEN
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
   90 CONTINUE
      RETURN
      END
