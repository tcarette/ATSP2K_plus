*
*     ------------------------------------------------------------------
*       C F G N 1
*     ------------------------------------------------------------------
*
      SUBROUTINE CFGN1(NCLOSD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NWD=94,NWCD=20)
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
      DATA iaster,iblnk/4h*   ,4h    /
*
      CALL CFGO1 (NCFG,MAXORB,NCLOSD,QIAJCMP,QLJCOMP,QNJCOMP,
     :                 QNOC,QNELCSH,QNOCORB,QJ1,QIAJCLD,QLJCLSD)
      NWF = MAXORB
 
*
*     In this version, all orbitals are assumed orthogonal!
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
      RETURN
      END
