!
!     ---------------------------------------------------------------
!        I N I T R
!     ---------------------------------------------------------------
!
!
      SUBROUTINE INITR 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE PARAM_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NOD = 220 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
! *****  SET THE COMMONLY USED DOUBLE PRECISION CONSTANTS
!
      D0 = 0.D0 
      D1 = 1.D0 
      D2 = 2.D0 
      D3 = 3.D0 
      D4 = 4.D0 
      D5 = 1.D0/2.D0 
      D6 = 6.D0 
      D8 = 8.D0 
      D10 = 10.D0 
      D12 = 12.D0 
      D16 = 16.D0 
      D30 = 30.D0 
!
! *****  SET THE STARTING POINT, STEP SIZE, AND RELATED PARAMETERS
!
      RHO = -4.D0 
      H = 1./16.D0 
      H1 = H/1.5 
      H3 = H/3. 
      CH = H*H/12. 
      EH = DEXP((-H)) 
      NO = NOD 
      ND = NO - 2 
!
! *****  SET THE FINE-STRUCTURE CONSTANT
!
      FINE = 0.25D0/137.036D0**2 
      RETURN  
      END SUBROUTINE INITR 
