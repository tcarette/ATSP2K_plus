*
*     -------------------------------------------------------------
*       C F G T S T
*     -------------------------------------------------------------
*
*     Modified by Gediminas Gaigalas,                September 1997
*
*
      SUBROUTINE CFGTST(NCFG,QLJCOMP,QNOC,QNELCSH,QNOCORB,QJ1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
*     THIS SUBROUTINE CHECKS ALL THE CONFIGURATION SET TO ENSURE THAT
*     IT SATISFIES ALL THE FOLLOWING CONDITIONS:
*        (1)  EACH CONFIGURATION HAS THE SAME NUMBER OF ELECTRONS
*        (2)  NO SUBSHELL HAS TOO MANY (.GT.2*(2*L+1))  ELECTRONS
*        (3)  THE ELECTRONS IN ANY ONE SUBSHELL ARE COUPLED TO FORM AN
*             ALLOWED TRIAD OF QUANTUM NUMBERS
*        (4)  THE TRIADS COUPLE TOGETHER IN AN ALLOWED WAY
*
*     IN THE EVENT OF AN ERROR, THE PROGRAM HALTS AT THE COMPLETION
*     OF THE CHECKING.  ANY NUMBER OF S, P, D  ELECTRONS ARE ALLOWED,
*     (BUT .LE.2*(2*L+1)), BUT ONLY UP TO TWO ELECTRONS, L >=3.
*     WHEN L>4, THE ONLY ALLOWED TERMS ARE THOSE FOR L=4.
*     A FILLED F-SHELL IS ALSO ALLOWED AS WELL AS A SINGLE ELECTRON
*     WITH L.GT.4
*
*
*     IMPLICIT INTEGER (Q)
      POINTER(QLJCOMP,LJCOMP(1))
      POINTER(QNOC,NOCCSH(1))
      POINTER(QNELCSH,NELCSH(8,1))
      POINTER(QNOCORB,NOCORB(8,1))
      POINTER(QJ1,J1QNRD(15,1))
*
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,JSC0,JSC(4)
      COMMON/TERMS/NROWS,ITAB(24),JTAB(24),NTAB(333)
CGGf
      COMMON /MT15/ M1(1),M2(7),M3(17),M4(47),M5(73)
      COMMON /MT67/ M7(119),M6(119)
      EXTERNAL TRMF,TRMF15,INITT
      DIMENSION IGMAX(14)
      DATA IGMAX/0,0,17,47,73,119,119,119,73,47,17,7,1,0/
CGGf
*
    5 FORMAT(/38H THE TRIAD OF QUANTUM NUMBERS OF SHELL,I3,17H IN CONFIG
     :URATION,I8,24H IS NOT A RECOGNIZED SET)
    7 FORMAT(/22H THE COUPLING OF SHELL,I3,17H IN CONFIGURATION,I3,
     : 38H RESULTS IN AN ILLEGAL COUPLING SCHEME)
   12 FORMAT(//41H CONFIGURATION DATA WRONG, PROGRAM HALTED//)
   15 FORMAT(/17H IN CONFIGURATION,I3,7H, SHELL,I3,28H CONTAINS TOO MANY
     : ELECTRONS)
   17 FORMAT(/14H CONFIGURATION,I3,68H INCLUDES A SHELL OF ANGULAR MOMEN
     :TUM L.GE.3 WITH TOO MANY ELECTRONS)
   18 FORMAT(/14H CONFIGURATION,I3,28H HAS AN INCORRECT NUMBER OF ,
     :        9HELECTRONS)
*
      IALLOW=1
      DO 1 I=1,NCFG
         NELSUM = 0
         N=NOCCSH(I)
         DO 2 J=1,N
            NA=NOCORB(J,I)
            NC=NELCSH(J,I)
            LQU=LJCOMP(NA)
            NELSUM = NELSUM + NC
            JD = J1QNRD(J,I)
            JA = MOD(JD,64)
            JD = JD/64
            JB = MOD(JD,64)
            JC = JD/64
            LQUMAX = 4*LQU + 2
            IF (NC .GT. LQUMAX) THEN
               WRITE(IWRITE,15) I,J
               IALLOW = 0
               GO TO 2
CGGf
            ELSE IF (LQU.EQ.3 .AND. NC.GT.2 .AND. NC.LT.14) THEN
            IMAXGG=IGMAX(NC)
            JCB=((JC-1)*100)+JB-1
              DO 14 IA =1,IMAXGG
                 IF (NC.EQ.3) THEN
                   LASTEL=M3(IA)
                 ELSEIF (NC.EQ.4) THEN
                   LASTEL=M4(IA)
                 ELSEIF (NC.EQ.5) THEN
                   LASTEL=M5(IA)
                 ELSEIF (NC.EQ.6) THEN
                   LASTEL=M6(IA)
                 ELSEIF (NC.EQ.7) THEN
                   LASTEL=M7(IA)
                 ELSEIF (NC.EQ.8) THEN
                   LASTEL=M6(IA)
                 ELSEIF (NC.EQ.9) THEN
                   LASTEL=M5(IA)
                 ELSEIF (NC.EQ.10) THEN
                   LASTEL=M4(IA)
                 ELSEIF (NC.EQ.11) THEN
                   LASTEL=M3(IA)
                 ELSEIF (NC.EQ.12) THEN
                   LASTEL=M2(IA)
                 ELSEIF (NC.EQ.13) THEN
                   LASTEL=M1(IA)
                 ELSE
                   STOP
                 ENDIF
                 LS=JTHN(LASTEL,1,10000)
                  IF(JCB.EQ.LS)THEN
                    NR=JTHN(LASTEL,4,100)
                    IF(JA.EQ.NR)GO TO 21
                  ENDIF
   14         CONTINUE
              GO TO 24
            ELSE IF (LQU.GT.4.AND.NC.GT.2) THEN
               WRITE(IWRITE,17) I
               IALLOW = 0
               GO TO 2
CGGf
            ELSE IF (NC .EQ. 1) THEN
               IF (JA.EQ.1 .AND. JB.EQ.(2*LQU+1) .AND. JC.EQ.2) GO TO 21
            ELSE
               IF (LQU .GT. 4 .AND. NC .EQ. 2) LQU = 4
               IF (NC .EQ. LQUMAX) THEN
                  NROW = 2
CGGf
               ELSEIF (NC .EQ. 0) THEN
                  NROW = 2
CGGf
               ELSE
                  NROW = NTAB1(NC+1,LQU+1)
               END IF
               I1 = ITAB(NROW)
               I2 = JTAB(NROW)
               DO 4 IA = 1,I1
                  I3 = I2+3*IA-1
                  IF (JB .EQ. NTAB(I3)) THEN
                     I3 = I3+1
                     IF (JC .EQ. NTAB(I3)) THEN
                        I3 = I3-2
                        IF (JA .EQ. NTAB(I3)) GO TO 21
                     END IF
                  END IF
    4          CONTINUE
            END IF
CGGf
   24       IALLOW = 0
CGGf
            WRITE(IWRITE,5) J,I
            GO TO 2
*
*     CHECK ON THE COUPLING OF THE TRIADS
*
   21       IF (N.GT.1 .AND. J.GT.1) THEN
               J2 = N+J-1
               J1 = J2-1
               IF (J.EQ.2) J1 = 1
               JE = J1QNRD(J1,I)/64
               JD = MOD(JE,64)
               JE = JE/64
               JG = J1QNRD(J2,I)/64
               JF = MOD(JG,64)
               JG = JG/64
               IF (JF.GE.(JB+JD) .OR. JF.LE.IABS(JB-JD) .OR.
     :             JG.GE.(JC+JE) .OR. JG.LE.IABS(JC-JE) .OR.
     :             MOD(JC+JE-JG,2).EQ.0 ) THEN
                   WRITE(IWRITE,7) J,I
                   IALLOW = 0
               END IF
            END IF
    2    CONTINUE
         IF (I .EQ. 1) THEN
            NELCS = NELSUM
         ELSE IF (NELSUM .NE. NELCS) THEN
            WRITE(IWRITE,18) I
            IALLOW = 0
         END IF
    1 CONTINUE
      IF (IALLOW .EQ. 0) THEN
         WRITE(IWRITE,12)
         STOP
      END IF
      END
