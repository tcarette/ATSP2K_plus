!
!     -------------------------------------------------------------
!      C A L C U L A T I O N
!     -------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!
      LOGICAL FUNCTION CALCULATION (KL) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE OPERAT_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:48:27  11/16/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: KL 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      CALCULATION = .FALSE. 
      IF (ICOLOM == 1) THEN 
         CALCULATION = .TRUE. 
      ELSE IF (IORBORB == 1) THEN 
         CALCULATION = .TRUE. 
      ELSE IF (ISOTOP == 1) THEN 
         IF (KL == 1) CALCULATION = .TRUE. 
      ENDIF 
      RETURN  
      END FUNCTION CALCULATION 
