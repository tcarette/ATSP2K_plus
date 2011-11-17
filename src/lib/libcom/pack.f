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
