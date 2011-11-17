*     ------------------------------------------------------------------
*
*       COMP -- A PROGRAM FOR PRINTING DOMINANT CONTRIBUTIONS
*
*       by C. Froese Fischer
*          Vanderbilt University
*          Nashville, TN 37235 USA
*
*       July, 1984
*
*-----------------------------------------------------------------------
*
      PROGRAM COMP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
ctc 14/10/2009
*      PARAMETER(NOD=220,NWD=60,NWD2=60,NCD=40000,NCD2=5*NCD)
      PARAMETER(NOD=220,NWD=60,NWD2=60,NCD=400000,NCD2=5*NCD)
ctc
      INTEGER Q(5),JV(NCD2),JL(NCD2),IP(NCD2)
      CHARACTER*1 ASTER,CH1,CH2
      CHARACTER*3 FIN1
      CHARACTER*3 EL(NWD),COUPLE(15),FINT(NCD2),elc(8)
      CHARACTER*6 ATOM,TERM,ALABEL
      CHARACTER*24 CFILE, NAME
      CHARACTER*64  CONFIG(NCD2), LAB1, LINE
      COMMON /STATE/LENGTH(NCD2),NCFG(2),NCF(2),IST(2)
      DIMENSION ET(100),COEF(NCD2,100)
      PARAMETER (IREAD=5,IWRITE=6)
      DATA ASTER/'*'/
 
*
      IU = 7
      DO 100 ICASE = 1,1
*
*  *****  DETERMINE INFORMATION ABOUT FILES
*
      iarg = iargc()
      if ( iarg .eq. 0) then
        WRITE(0,*) 'Enter name of .c file'
        READ(IREAD,'(A)') NAME
      else
        call getarg(1,NAME)
      end if
      jend = index(NAME,' ')
      CFILE = NAME(1:JEND-1)//'.c'
      OPEN(UNIT=7,FILE=CFILE,STATUS='OLD')
      M = 1
      WRITE(0,*) ' What tolerance? :'
      READ(IREAD,*) TOL
*
*  *****  READ THE CONFIGURATIONS
*
!     READ(IU,10) ATOM,ALABEL,ET(1),NCLOSD,NWF,NCFG(M)
!      print*, ATOM,ALABEL,ET(1),NCLOSD,NWF,NCFG(M)
!  10 FORMAT(3X,2A6,F14.7,1X,I3,I4,I5)
!     Read(iu,'(20(1x,A3):)') (el(i),i=1,nclosd)
!     IF (nwf .gt. nclosd)
!    :   read(iu,'(20(1X,A3):)') (el(i),i=nclosd+1,nwf)

!     The format has changed.  Read two lines:
      Read (IU, '(3X,2A6)') Atom, Term
      Read (IU, '(A)') Alabel

   15 READ(IU,'(8(1x,A3,1x,I2,1x),f10.8)',END=18) 
     :     (ELc(K),Q(K),K=1,8),COEF(M,1)
      READ(IU,'(15(1X,A3))',END=18) (COUPLE(J),J=1,15)
      NOCC = 0
16    IF (ELc(NOCC+1) .NE. '   ' ) THEN
         NOCC = NOCC+1
         IF (NOCC .LT. (8)) GO TO 16
      END IF
      IF (NOCC .EQ. 0) GO TO 18
      CALL PACK(NOCC,ELc,Q,COUPLE,LAB1)
 
*
*       1. Separate the final term
*
            J = INDEX(LAB1,' ')-1
           CH1 = LAB1(J:J)
           IF (CH1 .GE. '0' .AND. CH1 .LE. '9') THEN
                J = J-4
            ELSE
                J = J-3
            ENDIF
            FIN1 = LAB1(J+2:J+3)
            LAB1(J+1:J+4) = '    '
 
*
*       2. Delete set subscriP2
*
*           CH1 = LAB1(J:J)
*           IF (CH1.GE.'0' .AND. CH1.LE.'9') LAB1(J:J) = ' '
 
*
*       3. If after removing the final term, there are no other
*   intermediate couplings prefaced by '_' and the last coupling
*   is the same as the final term, then the coupling for the final
*   term is omitted .
*
            IF (INDEX(LAB1,'_') .EQ. 0) THEN
                CH2 = LAB1(J-2:J-1)
                IF (CH2 .EQ. FIN1) LAB1(J-2:J-1) = '  '
            ENDIF
      CONFIG(M) = LAB1
      K = 64
   17 IF (CONFIG(M)(K:K) .EQ. ' ') THEN
         K = K-1
         GO TO 17
      END IF
      LENGTH(M) = K
      FINT(M) = FIN1
      M = M+1
      IF (M.GT.(NCD2)) THEN
         WRITE(0,'(A,I4)')
     :       ' TOO MANY CONFIGURATIONS:  MAXIMUM IS ', (NCD2)-1
      END IF
      GO TO 15
   18 NCF(ICASE)= M-1
      CLOSE(UNIT=7)
100   CONTINUE
 
      WRITE(0,'(A/5X,A/5X,A/5X,A/A)') ' Compositions from:',
     :     '1  name.c','2  name.l','3  name.j',' Enter selection'
      READ(IREAD,*) ICASE
      IF (ICASE .EQ. 1) THEN
	 NCFG(1) = NCF(1)
	 IP(1) = 1
	 JL(1) = 1
	 IST(1) = 1
         CALL OUTPUT(1,O,JL,CONFIG,FINT,COEF,IP,TOL,ET)
      ELSE
      IF (ICASE .EQ. 2) THEN
	 CFILE = NAME(1:JEND-1)//'.l'
      ELSE 
	 CFILE = NAME(1:JEND-1)//'.j'
      END IF
*
*   ****  If LS format, read sets of coefficients
*
      IS = 1
      DO 200 ICASE = 1,1
      IL = NCF(ICASE)
      OPEN(UNIT=7,FILE=CFILE,STATUS='OLD')
      READ(7,'(A6,7X,F6.1,7X,I4,9X,I7)') ATOM,Z,NEL,NCFG(ICASE)
65    READ(7,'(//8X,I4,10X,I4)',END=190) JVV,MFOUND
      DO 63 III = 1,MFOUND
        READ(7,64) JL(IS),ET(IS),(COEF(I,IS),I=1,NCFG(ICASE))
64      FORMAT(/I6,F16.8/(7F11.8))
        IF (ICASE .EQ. 2) JL(IS) = JL(IS) + NCF(1)
        JV(IS) = JVV
        IS = IS+1
	if (is .gt. 100) then
	   write(0,*) ' Current dimension exeeded: Max=100'
	   stop
	end if
63    CONTINUE
      GO TO 65
190   IST(ICASE) = IS-1
200   CONTINUE
      IS = IS-1
      CALL SORT(IS,ET,IP)
      CALL OUTPUT(IS,JV,JL,CONFIG,FINT,COEF,IP,TOL,ET)
      END IF
      END
*
*     ------------------------------------------------------------------
*       OUTPUT
*     ------------------------------------------------------------------
*
      SUBROUTINE OUTPUT(IS,JV,JL,CONFIG,FINT,COEF,IP,TOL,ET)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
ctc 14/10/2009
*      PARAMETER(NOD=220,NWD=60,NWD2=60,NCD=40000,NCD2=5*NCD)
      PARAMETER(NOD=220,NWD=60,NWD2=60,NCD=400000,NCD2=5*NCD)
ctc
      PARAMETER (IREAD=5,IWRITE=6)
      DIMENSION ET(*), COEF(NCD2,100)
      COMMON /STATE/LENGTH(NCD2),NCFG(2),NCF(2),IST(2)
      CHARACTER*64 CONFIG(*)
      CHARACTER*3 FINT(*), A*1
      INTEGER INDEX(NCD2),JV(NCD2),JL(NCD2),IP(*)

      LOGICAL FIRST,SECOND
      DATA A/'&'/
! Configuration  Term  ****     **************
! 2s(2).2p(6).3s(2)  1S     .0     -1621.08109002

      DO 100  II = 1,IS
         JS = IP(II)
         JSL = JL(JS)
         ICASE = 1
         IF (JSL .GT. NCF(1)) ICASE = 2
         K = LENGTH(JSL)
      write(iwrite,'(//1X,A30,2X,A4,2X,A4X,A15)')
     :                 'Configuration', 'Term', 'J','Total Energy'
           WRITE(IWRITE,'(1X,A30,2X,A4,2X,F4.1,5X,F14.8)')
     :            CONFIG(JSL)(1:K),FINT(JSL),JV(JS)/2.,ET(II)
      DO 1 I = 1,NCFG(ICASE)
         INDEX(I) = I
1     CONTINUE
      DO 2 I = 1,NCFG(ICASE)
         JP = I
         DO 3 J = I+1,NCFG(ICASE)
            IF (ABS(COEF(J,JS)) .GT. ABS(COEF(JP,JS))) JP = J
3        CONTINUE
         TEMP = COEF(I,JS)
         COEF(I,JS) = COEF(JP,JS)
         COEF(JP,JS) = TEMP
         ITEMP = INDEX(I)
         INDEX(I) = INDEX(JP)
         INDEX(JP) = ITEMP
2     CONTINUE
        IB = 0
        IF (JS .GT. IST(1)) IB = NCF(1)
          FIRST = .TRUE.
          SECOND = .FALSE.
      DO 10 I = 1,NCFG(ICASE)
         J = INDEX(I) + IB
         K = LENGTH(J)
         IF (ABS(COEF(I,JS)) .GT. TOL) THEN
            IF (FIRST) THEN
                   WRITE(IWRITE,'(1X,F12.7,2X,A50,2X,A)')
     :                  COEF(I,JS),CONFIG(J)(1:K),FINT(J)
                FIRST = .FALSE.
                SECOND = .TRUE.
            ELSE IF (SECOND) THEN
                   WRITE(IWRITE,'(1X,F12.7,2X,A50,2X,A)')
     :                  COEF(I,JS),CONFIG(J)(1:K),FINT(J)
                SECOND = .FALSE.
            ELSE
                   WRITE(IWRITE,'(1X,F12.7,2X,A50,2X,A)')
     :                  COEF(I,JS),CONFIG(J)(1:K),FINT(J)
               SECOND = .TRUE.
            END IF
 
        END IF
10      CONTINUE
100    CONTINUE
      END
*
*     ------------------------------------------------------------------
*       SORT
*     ------------------------------------------------------------------
*
        SUBROUTINE SORT(M,ET,IP)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION ET(*),IP(*)
 
        DO 1 I = 1,M
           IP(I) = I
1       CONTINUE
 
        DO 10 I = 1,M-1
           JP = I
           DO 12 J = I+1,M
              IF (ET(J) .LT. ET(JP)) JP = J
12         CONTINUE
           TEMP = ET(I)
           ET(I) = ET(JP)
           ET(JP) = TEMP
           ITEMP = IP(I)
           IP(I) = IP(JP)
           IP(JP) = ITEMP
10      CONTINUE
        END
 

*     -------------------------------------------------------------
*       L V A L
*     -------------------------------------------------------------
*
*     Modified by Gediminas Gaigalas,                September 1997
*
*
      INTEGER FUNCTION LVAL(SYMBOL)
      CHARACTER*1 SYMBOL
      CHARACTER*26 SET
      DATA SET/'spdfghiklmnoqSPDFGHIKLMNOQ'/
*
      LOCATE = INDEX(SET,SYMBOL)
      IF ( LOCATE .LE. 13) THEN
            LVAL = LOCATE - 1
         ELSE
            LVAL = LOCATE - 14
      ENDIF
      RETURN
      END
*
*     ------------------------------------------------------------------
*       P A C K
*     ------------------------------------------------------------------
*  Subroutine written by Bin LIU
*
*  Rules for encoding
*  1. All blanks deleted
*  2. If Qi=1, omit Qi
*  3. If Qi=1 or Qi>=4l+1, omit ALFAi
*  4. If i=1 or (Qi=4l+2 and i<>m), insert '.'; else _BETAi.
*
      SUBROUTINE PACK (M, EL, Q, COUPLE, STR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 
          INTEGER FULL,Q(*),CONST
          CHARACTER*3 EL(*),COUPLE(*)
          CHARACTER   STR*64,CH1*1,CH3*3
*
*  FULL   :  4l+2
*  CONST  :   constant for converting lowercase to uppercase
*  CH*    :  temporary variables
*
          CONST = ICHAR('a') - ICHAR('A')
          STR=' '
          J = 1
*
*   -----  begin to encode  -----
*
          DO 100 I=1,M
            IF (Q(I) .EQ. 0) GO TO 100
            IF (EL(I)(1:1) .EQ. ' ') THEN
                STR(J:J+1) = EL(I)(2:3)
                K = 2
            ELSE
                STR(J:J+2) = EL(I)
                K = 3
                IF (EL(I)(3:3) .EQ. ' ') K = 2
            END IF
            CH1 = STR(J+1:J+1)
            IF ((CH1.GE.'A') .AND. (CH1.LE.'Z'))
     :            STR(J+1:J+1)=CHAR(ICHAR(CH1)+CONST)
            FULL=4*LVAL(CH1)+2
*
*  -----  convert Qi into character  -----
*
            J = J+K
            N=Q(I)
            IF (N .GT. 14) STOP 'Too many electrons'
*
*  -----  If Qi<>1, add Qi
*           If Qi<4l+1, add TERMi for the shell -----
*
            IF (N .NE. 1) THEN
               STR(J:J) = '('
               J = J + 1
               IF (N .GT. 9) THEN
                    STR(J:J) = '1'
                    J = J + 1
                    N = N - 10
               END IF
               STR(J:J+1) = CHAR(ICHAR('0')+N)//')'
               J = J + 2
               IF (N .LT. FULL-1 .AND. M .NE. 1) THEN
                   CH3=COUPLE(I)
                   CH1= CH3(2:2)
                   IF (CH1.GE.'a' .AND. CH1.LE.'z')
     :                 CH3(2:2) = CHAR(ICHAR(CH1)-CONST)
                   STR(J:J+2) = CH3
                   J = J+3
               ENDIF
            ENDIF
*
*  -----  If i=1 or Qi=4l+2 and i<>m,
*           insert '.'; else _RESULTANTi.  -----
*
            IF ((I .NE. 1 .AND. N .NE. FULL) .OR. M .EQ. I) THEN
                CH3=COUPLE(M+I-1)
                CH1 = CH3(2:2)
                IF (CH1.GE.'a' .AND. CH1.LE.'z')
     :                  CH3(2:2) = CHAR(ICHAR(CH1)-CONST)
                K = 2
                IF (M .EQ. 1) K = 3
                STR(J:J+K) = '_'//CH3(1:K)
                J = J + K + 1
            ENDIF
            IF (I .NE. M) THEN
                STR(J:J)='.'
                J = J + 1
            ENDIF
100       CONTINUE
*
*>>>>>    Because of a compiler error on the SUN, the following is
*         needed to have the string printed correctly.
*         STR = STR(1:J)
          RETURN
        END
