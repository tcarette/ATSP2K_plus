*    --------------------------------------------------------------
*            R E F O R M
*    --------------------------------------------------------------
*
*
      SUBROUTINE REFORM(STR1,STR2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*40 STR1,STR2,BLANK
      DATA BLANK/'   '/
*
    1 I = 0
      STR2 = BLANK
      IS = 0
    2 JS = INDEX(STR1(IS+1:),'(')
      IF (JS .NE. 0) THEN
         IF (JS .GT. 5) GO TO 10
         I = I+5
         STR2(I-JS+1:I) = STR1(IS+1:IS+JS)
         IS = IS + JS
         JS = INDEX(STR1(IS+1:),')')
         IF (JS .EQ. 0 .OR. JS .GT. 3) GO TO 10
         I = I+3
         STR2(I-JS+1:I) = STR1(IS+1:IS+JS)
         IS = IS + JS
         GO TO 2
      END IF
      RETURN
   10 PRINT *,' Error in ',STR1,': Re-enter'
      READ '(A)',STR1
      GO TO 1
      END
