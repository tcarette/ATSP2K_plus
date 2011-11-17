!    --------------------------------------------------------------
!            R E F O R M
!    --------------------------------------------------------------
!
!
      SUBROUTINE REFORM(STR1, STR2) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER  :: STR1*40 
      CHARACTER , INTENT(OUT) :: STR2*40 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IS, JS 
      CHARACTER :: BLANK*40 
!-----------------------------------------------
      DATA BLANK/ '   '/  
!
    1 CONTINUE 
      I = 0 
      STR2 = BLANK 
      IS = 0 
    2 CONTINUE 
      JS = INDEX(STR1(IS+1:),'(') 
      IF (JS /= 0) THEN 
         IF (JS > 5) GO TO 10 
         I = I + 5 
         STR2(I-JS+1:I) = STR1(IS+1:IS+JS) 
         IS = IS + JS 
         JS = INDEX(STR1(IS+1:),')') 
         IF (JS==0 .OR. JS>3) GO TO 10 
         I = I + 3 
         STR2(I-JS+1:I) = STR1(IS+1:IS+JS) 
         IS = IS + JS 
         GO TO 2 
      ENDIF 
      RETURN  
   10 CONTINUE 
      WRITE (6, *) ' Error in ', STR1, ': Re-enter' 
      READ '(A)' , STR1 
      GO TO 1 
      END SUBROUTINE REFORM 
