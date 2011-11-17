*
*     ------------------------------------------------------------------
*	C F G T S T
*     ------------------------------------------------------------------
*
      SUBROUTINE CFGTST(NCFG,LJCOMP,NOCCSH,NELCSH,NOCORB,J1QNRD,NCD)
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
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,JSC0,JSC(3)
      COMMON/TERMS/NROWS,ITAB(24),JTAB(24),NTAB(333)
*
      DIMENSION LJCOMP(*),NOCCSH(NCD),NELCSH(5,NCD),NOCORB(5,NCD),
     :          J1QNRD(9,NCD)
*
    5 FORMAT(/38H THE TRIAD OF QUANTUM NUMBERS OF SHELL,I3,17H IN CONFIG
     :URATION,I3,24H IS NOT A RECOGNIZED SET)
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
            LQU=LJCOMP(NA)
            NC=NELCSH(J,I)
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
            ELSE IF ((LQU.EQ.3 .AND. NC.GT.2 .AND. NC.LT.14) .OR.
     :          (LQU.GT.4.AND.NC.GT.2)) THEN
               WRITE(IWRITE,17) I
               IALLOW = 0
               GO TO 2
            ELSE IF (NC .EQ. 1) THEN
               IF (JA.EQ.1 .AND. JB.EQ.(2*LQU+1) .AND. JC.EQ.2) GO TO 21
            ELSE
               IF (LQU .GT. 4 .AND. NC .EQ. 2) LQU = 4
               IF (NC .EQ. LQUMAX) THEN
                  NROW = 2
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
            IALLOW = 0
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
