!
!     ------------------------------------------------------------------
!     F I N D
!     ------------------------------------------------------------------
!
      CHARACTER*(3) FUNCTION FIND (I, OF, EL) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:26:11  11/16/01  
!...Switches:                     
!
!  ---  THIS ROUTINE FINDS ELECTRONS IN ONE OF THREE LISTS
!
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NWD = 128 
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I 
      CHARACTER , INTENT(IN) :: OF(NWD,2)*3 
      CHARACTER , INTENT(IN) :: EL(NWD)*3 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
      IF (I <= NWD) THEN 
         FIND = EL(I) 
      ELSE IF (I <= NWD*2) THEN 
         FIND = OF(I-NWD,1) 
      ELSE 
         FIND = OF(I-2*NWD,2) 
      ENDIF 
      RETURN  
      END FUNCTION FIND 
