!
!     ------------------------------------------------------------------
!             E P T R
!     ------------------------------------------------------------------
!
!       Determines the position of the electron in the electron list
!
      SUBROUTINE EPTR(EL, ELSYMB, IEL, J2) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE PARAM_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(OUT) :: IEL 
      CHARACTER , INTENT(IN) :: ELSYMB*3 
      CHARACTER , INTENT(IN) :: EL(*)*3 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: IWRITE = 6 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J2, I, J1 
      CHARACTER :: BL*3 
!-----------------------------------------------
      DATA BL/ '   '/  
      J2 = 0 
!
! ***** SEARCH ELECTRON LIST FOR LSYMB
!
      IF (ELSYMB == BL) THEN 
         IEL = 0 
         RETURN  
      ENDIF 
      DO I = 1, NWF 
         IF (EL(I) /= ELSYMB) CYCLE  
         IEL = I 
         RETURN  
      END DO 
      IEL = -1 
      WRITE (IWRITE, 20) ELSYMB 
   20 FORMAT(/,10X,A3,' NOT FOUND IN ELECTRON LIST') 
      J2 = 1 
      RETURN  
      END SUBROUTINE EPTR 
