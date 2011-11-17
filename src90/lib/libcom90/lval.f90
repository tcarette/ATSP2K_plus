!     -------------------------------------------------------------
!       L V A L
!     -------------------------------------------------------------
!
!     Modified by Gediminas Gaigalas,                September 1997
!
!
      INTEGER FUNCTION LVAL (SYMBOL) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER , INTENT(IN) :: SYMBOL 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LOCATE 
      CHARACTER :: SET*26 
!-----------------------------------------------
      DATA SET/ 'spdfghiklmnoqSPDFGHIKLMNOQ'/  
!
      LOCATE = INDEX(SET,SYMBOL) 
      IF (LOCATE <= 13) THEN 
         LVAL = LOCATE - 1 
      ELSE 
         LVAL = LOCATE - 14 
      ENDIF 
      RETURN  
      END FUNCTION LVAL 
