!
!     -------------------------------------------------------------
!       N U M V A L
!     -------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      INTEGER FUNCTION NUMVAL (SYMBOL) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  20:44:41  11/14/01  
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
      CHARACTER :: SET*13 
!-----------------------------------------------
      DATA SET/ '0123456789 aA'/  
!
      LOCATE = INDEX(SET,SYMBOL) 
      IF (LOCATE <= 10) THEN 
         NUMVAL = LOCATE - 1 
      ELSE IF (LOCATE == 11) THEN 
         NUMVAL = 0 
      ELSE 
         NUMVAL = 10 
      ENDIF 
      RETURN  
      END FUNCTION NUMVAL 
