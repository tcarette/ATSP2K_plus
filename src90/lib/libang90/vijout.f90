!
!     ------------------------------------------------------------------
!       V I J O U T
!     ------------------------------------------------------------------
!
      SUBROUTINE VIJOUT(JA, JB) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE INFORM_C 
      USE DEBUG_C 
      USE MEDEFN_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:26:42  11/16/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: JA 
      INTEGER , INTENT(IN) :: JB 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I2HSH, I, J, K 
!-----------------------------------------------
!
!     THIS SUBROUTINE IS ENTERED ONLY IF IBUG2 IS GREATER THAN ZERO
!
! --- PRINT OUT OF QUANTUM NUMBERS AND COUPLING SCHEMES FOR EACH
!     MATRIX ELEMENT AS DEFINED BY SETUP
!
    5 FORMAT(/,/,' L.H.S. OF HAMILTONIAN MATRIX ELEMENT DEFINED BY') 
    6 FORMAT(/,/,' R.H.S. OF HAMILTONIAN MATRIX ELEMENT DEFINED BY') 
    7 FORMAT('1(CONFIG ',I2,'/V/CONFIG ',I2,')') 
    8 FORMAT(/,' NJ,LJ ',10(I6,I3)) 
    9 FORMAT(/,' NOSH ',10I4) 
   10 FORMAT(' J1QN ',10(I5,2I3)) 
      I2HSH = 2*IHSH - 1 
      WRITE (IWRITE, 7) JA, JB 
      WRITE (IWRITE, 8) (NJ(I),LJ(I),I=1,IHSH) 
      WRITE (IWRITE, 5) 
      WRITE (IWRITE, 9) (NOSH1(J),J=1,IHSH) 
      WRITE (IWRITE, 10) ((J1QN1(J,K),K=1,3),J=1,I2HSH) 
      WRITE (IWRITE, 6) 
      WRITE (IWRITE, 9) (NOSH2(J),J=1,IHSH) 
      WRITE (IWRITE, 10) ((J1QN2(J,K),K=1,3),J=1,I2HSH) 
      RETURN  
      END SUBROUTINE VIJOUT 
