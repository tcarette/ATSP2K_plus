*
*     ------------------------------------------------------------------
*	C F G I N 2
*     ------------------------------------------------------------------
*
      SUBROUTINE CFGIN2(MCFG,KCFG,LORTH,INPUT)
*
*	Read two sets of configurations and determine the orthogonality
*       conditions between them
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NCD=2000,NCD4=4*NCD)
      CHARACTER*3 EL(2*NWD), ELC(NWD), JAJCMP(2*NWD,3)*1
      CHARACTER INPUT(2)*24,HEADI*72,HEADF*72,HEADER*72
      LOGICAL LORTH
      COMMON/INFORM/ IREADI,IWRITE,IOUT,IREADF,ISC(7)
      COMMON/STATES/NCFG,MAXORB,IAJCMP(2*NWD),LJCOMP(2*NWD),
     :NJCOMP(2*NWD),NOCCSH(NCD4),NELCSH(5,NCD4),NOCORB(5,NCD4),
     :J1QNRD(9,NCD4)
      COMMON/NOR/NCOM,NORBI,NORBF,IWAR
*
    3 FORMAT(18(1X,A3))
    7 FORMAT(A72)
   22 FORMAT(// 7H STATE ,' (WITH',I3,' CONFIGURATIONS):'/1H ,31(1H-)/)
   23 FORMAT(/10H THERE ARE,I3,21H ORBITALS AS FOLLOWS://
     1 5X,21(1X,A3):/5X,21(1X,A3))
*
*     The "readonly" option is needed on some computers when the
*     two files are in fact the same.  Others ignore the option 
*     which is OK in most cases.
*     Microsoft Fortran requires "READ" instead of "READONLY"
*     OPEN(UNIT=1,FILE=INPUT(1),STATUS='OLD',READONLY)
*     OPEN(UNIT=2,FILE=INPUT(2),STATUS='OLD',READONLY)
      OPEN(UNIT=1,FILE=INPUT(1),STATUS='OLD')
      OPEN(UNIT=2,FILE=INPUT(2),STATUS='OLD')
*
* --- ANALYZE INITIAL AND FINAL STATE DATA
*
      CALL ANALY2(NCLOSI,NCLOSF,MCFG,KCFG,EL,LORTH)
      REWIND(UNIT=IREADI)
      REWIND(UNIT=IREADF)
*
      MAXORB = NCOM + NORBI + NORBF
*   SET UP THE ELECTRONS
*
      READ(EL,'(A3)') (IAJCMP(I),I=1,MAXORB)
      READ(EL,'(3A1)')((JAJCMP(I,J),J=1,3),I=1,MAXORB)
*
*   SET UP OF LJCOMP
*
      DO 60 I = 1,MAXORB
      IF (JAJCMP(I,1) .EQ. ' ') THEN
         JAJCMP(I,1) = JAJCMP(I,2)
         JAJCMP(I,2) = JAJCMP(I,3)
         JAJCMP(I,3) = ' '
      ENDIF
      LJCOMP(I) = LVAL(JAJCMP(I,2))
      NJCOMP(I) = ICHAR(JAJCMP(I,1)) - ICHAR('1') + 1
   60 CONTINUE
*
* ---- CHECK COMMON CLOSED SHELLS
*
      IF (NCLOSI .NE. NCLOSF)
     :   STOP ' Common closed shells not the same in the two states'
*
      READ(IREADI,7) HEADI
      READ(IREADF,7) HEADF
      HEADER = HEADI(1:34)//'=>'//HEADF(1:34)
      WRITE(IOUT,7) HEADER
*
* --- CHECK CLOSED SHELLS FURTHER
*
      READ(IREADI,3) (ELC(I),I=1,NCLOSI)
      READ(IREADF,3) (EL(I),I=1,NCLOSF)
      DO 1 I = 1,NCLOSF
         J = 1
    2    IF (EL(I) .NE. ELC(J) ) THEN
            J = J+1
            IF (J .LE. NCLOSI) THEN
               GO TO 2
              ELSE
               STOP ' Common closed sub-shells not the same'
            END IF
         END IF
    1 CONTINUE
*
*  MAXORB < 2*NWD+1... LINKED TO JAJCMP(2*NWD,3)
*                            IAJCMP(2*NWD)
*                            LJCOMP(2*NWD)
*                            NJCOMP(2*NWD)
*  THE DIMENSION OF IORTH IS GIVEN BY THE PRODUCT OF THE ALLOWED
*  NORBI AND NORBF, I.E. ACTUALLY NWD by NWD = NWD^2
*
*
*   GET INITIAL STATE CONFIGURATIONS
*
      CALL GSTATE(1,MCFG)
*
*
      MCFG1 = MCFG + 1
      NCFG = MCFG + KCFG
      CALL GSTATE(MCFG1,NCFG)
*
*  ---  CHECK THE DATA
*
      CALL CFGTST(NCFG,LJCOMP,NOCCSH,NELCSH,NOCORB,J1QNRD,NCD4)
      RETURN
      END
